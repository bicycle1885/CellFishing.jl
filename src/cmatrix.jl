# Compressed Matrix
struct CMatrix{T<:Integer} <: AbstractMatrix{T}
    data::Vector{UInt8}
    size::Tuple{Int,Int}
end

CMatrix(counts::Matrix; kwargs...) = CMatrix{eltype(counts)}(blosc_compress(counts; kwargs...), size(counts))

Base.size(M::CMatrix) = M.size
Base.IndexStyle(::Type{<:CMatrix}) = IndexLinear()
Base.getindex(M::CMatrix, i::Integer) = getitem(M.data, eltype(M), i)
Base.getindex(M::CMatrix, ::Colon, j::Integer) = getitems(M.data, eltype(M), (j - 1) * size(M, 1) + 1, size(M, 1))

const BLOSC_MAX_OVERHEAD = 16  # BLOSC_MIN_HEADER_LENGTH

function blosc_compress(
        data::Matrix;
        level::Integer=8,
        shuffle::Bool=true,
        compressor::String="blosclz",
        blocksize::Integer=0,
        nthreads::Integer=1,
    )
    srcsize = sizeof(data)
    src = pointer(data)
    dstsize = srcsize + BLOSC_MAX_OVERHEAD
    dst = Vector{UInt8}(undef, dstsize)
    ret = ccall(
        (:blosc_compress_ctx, libblosc),
        Cint,
        (Cint, Cint, Csize_t,
         Csize_t, Ptr{Cvoid}, Ptr{Cvoid},
         Csize_t, Cstring,
         Csize_t, Cint),
        level, shuffle, sizeof(eltype(data)),
        srcsize, src, dst,
        dstsize, compressor,
        blocksize, nthreads)
    @assert ret ≥ 0
    resize!(dst, ret)
    return dst
end

function getitem(data::Vector{UInt8}, ::Type{T}, i::Integer) where T
    src = pointer(data)
    dst = Ref{T}()
    ret = ccall(
        (:blosc_getitem, libblosc),
        Cint,
        (Ptr{Cvoid}, Cint, Cint, Ref{T}),
        src, i - 1, 1, dst)
    @assert ret == sizeof(T)
    return dst[]
end

function getitems(data::Vector{UInt8}, ::Type{T}, i::Integer, n::Integer) where T
    src = pointer(data)
    dst = Vector{T}(undef, n)
    ret = ccall(
        (:blosc_getitem, libblosc),
        Cint,
        (Ptr{Cvoid}, Cint, Cint, Ref{T}),
        src, i - 1, n, dst)
    @assert ret == sizeof(T) * n
    return dst
end

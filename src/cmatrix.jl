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
function Base.getindex(M::CMatrix, ::Colon, cols::AbstractVector{<:Integer})
    n = length(cols)
    A = Matrix{eltype(M)}(undef, size(M, 1), n)
    for j in 1:n
        A[:,j] .= M[:,cols[j]]
    end
    return A
end

const BLOSC_MAX_OVERHEAD = 16  # BLOSC_MIN_HEADER_LENGTH

function blosc_compress(
        data::Matrix{T};
        level::Integer=8,
        shuffle::Bool=true,
        compressor::String="lz4hc",
        blocksize::Integer=size(data, 1) * sizeof(T),
        nthreads::Integer=1,
        mmap::Bool=false,
    ) where {T<:Integer}
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
        level, shuffle, sizeof(T),
        srcsize, src, dst,
        dstsize, compressor,
        blocksize, nthreads)
    @assert ret ≥ 0
    if mmap
        dst_mmap = Mmap.mmap(Vector{UInt8}, (ret,), shared=false)
        copyto!(dst_mmap, 1, dst, 1, ret)
        return dst_mmap
    else
        resize!(dst, ret)
        return dst
    end
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

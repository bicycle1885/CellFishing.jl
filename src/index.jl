# Cell Index
# ==========

# Runtime Statistics
# ------------------

mutable struct RuntimeStats
    preproc_time::UInt64
    search_time::UInt64
    rank_time::UInt64
    timer::UInt64
end

RuntimeStats() = RuntimeStats(0, 0, 0, 0)

function reset!(s::RuntimeStats)
    s.preproc_time = 0
    s.search_time = 0
    s.rank_time = 0
    s.timer = 0
    return s
end

@inline function tic!(s::RuntimeStats)
    t = time_ns()
    ns = t - s.timer
    s.timer = t
    return ns
end


# Locality-Sensitive Hash
# -----------------------

struct LSHash{T<:BitVec}
    projection::Matrix{Float32}
    hammingindex::HammingIndexes.HammingIndex{T}
end

bitvectype(::LSHash{T}) where T = T

function findneighbors!(neighbors::Matrix{Int},
                        Z::Vector{T},
                        X::Matrix{Float32},
                        lshash::LSHash{T}) where T
    @assert length(Z) == size(neighbors, 2) == size(X, 2)
    sketch!(Z, X, lshash)
    # search for k-nearest neighbors
    k = size(neighbors, 1)
    nns = HammingIndexes.NearestNeighbors(k)
    for j in 1:length(Z)
        HammingIndexes.findknn!(nns, Z[j], lshash.hammingindex)
        for i in 1:k
            neighbors[i,j] = nns.indexes[i]
        end
    end
    return neighbors, Z
end

# Draw sketches of X.
function sketch!(Z::Vector{T}, X::Matrix{Float32}, lshash::LSHash{T}) where T
    return sketch!(Z, X, lshash.projection)
end

function sketch!(Z::Vector{T}, X::Matrix{Float32}, P::Matrix{Float32}) where {T<:BitVec}
    n = size(X, 2)
    @assert n == length(Z)
    PX = P * X
    @assert size(PX, 1) == bitsof(T)
    for j in 1:n
        z::T = 0
        i = (j - 1) * bitsof(T) + 1
        @inbounds for k in 0:bitsof(T)-1
            z |= T(PX[i+k] ≥ 0) << k
        end
        Z[j] = z
    end
    return Z
end

function generate_random_projections(K::Integer, D::Integer, superbit::Integer)
    @assert 1 ≤ superbit ≤ min(K, D)
    q, r = divrem(K, superbit)
    Ps = Matrix{Float32}[]
    for _ in 1:q
        # generate an orthogonal random matrix.
        M = randn(Float32, D, superbit)
        push!(Ps, Matrix(qr!(M).Q)')
    end
    if r > 0
        M = randn(Float32, D, r)
        push!(Ps, Matrix(qr!(M).Q)')
    end
    return vcat(Ps...)
end


# CellIndex
# ---------

struct CellIndex
    # preprocessor
    preproc::Preprocessor
    # locality-sensitive hashes
    lshashes::Vector{LSHash{T}} where T
    # any user-provided metadata (e.g. cellnames)
    metadata::Any
    # compressed count matrix (optional)
    counts::Union{CMatrix{Int32},Nothing}
    # run-time statistics
    rtstats::RuntimeStats
end

nbits(index::CellIndex) = bitsof(bitvectype(index.lshashes[1]))
ncells(index::CellIndex) = length(index.lshashes[1].hammingindex)

function Base.show(io::IO, index::CellIndex)
    print(io, summary(index), "(<#cells=$(ncells(index)), hash=$(nbits(index))×$(length(index.lshashes))>)")
end

"""
    CellIndex{T}(counts, features; <keyword arguments>...)

Create a cell index from a count matrix `counts` using `features`.

Arguments
---------

- `counts`: transcriptome expression matrix (features x cells) [required].
- `features`: features used to compare expression profiles [required].
- `scalefactor=1.0e4`: the scale factor of library sizes.
- `n_dims=50`: the number of dimensions after PCA.
- `transformer=:log1p`: the variance-stabilizing transformer (`:log1p` or `:ftt`).
- `randomize=true`: to use randomized SVD or not.
- `normalize=true`: to normalize library sizes or not.
- `standardize=true`: to standardize features or not.
- `metadata`: arbitrary metadata.
- `n_bits=128`: the number of bits (64, 128, 256, or 512).
- `n_lshashes=4`: the number of locality-sensitive hashes.
- `superbit=min(n_dims, n_bits)`: the depth of super-bits.
- `index=true`: to create bit index(es) or not.
- `keep_counts=false`: to keep filtered counts in the index object.
"""
function CellIndex(
        # required arguments
        counts::AbstractMatrix{<:Real},
        features::Features;
        # additional data
        featurenames=nothing,
        metadata=nothing,
        # parameters for preprocessing
        scalefactor::Real=1.0e4,
        n_dims::Integer=50,
        transformer::Symbol=:log1p,
        randomize::Bool=true,
        normalize::Bool=true,
        standardize::Bool=true,
        # parameters for LSH
        n_bits::Integer=128,
        n_lshashes::Integer=4,
        superbit::Integer=min(n_dims, n_bits),
        index::Bool=true,
        keep_counts::Bool=false,
       )
    # check arguments
    m, n = size(counts)
    if nfeatures(features) != m
        throw(ArgumentError("mismatching features size"))
    end
    if nselected(features) == 0
        throw(ArgumentError("no selected features"))
    end
    if !(1 ≤ n_dims ≤ m)
        throw(ArgumentError("invalid n_dims"))
    end
    if n_bits ∉ (64, 128, 256, 512)
        throw(ArgumentError("invalid n_bits"))
    end
    if !(1 ≤ superbit ≤ min(n_dims, n_bits))
        throw(ArgumentError("invalid superbit"))
    end
    if !(scalefactor > 0)
        throw(ArgumentError("invalid scalefactor"))
    end
    if !(1 ≤ superbit ≤ min(n_dims, n_bits))
        throw(ArgumentError("invalid superbit"))
    end
    if transformer ∉ (:log1p, :ftt)
        throw(ArgumentError("invalid transformer"))
    end
    # filter features
    featurenames = selectedfeatures(features)
    counts = Matrix(counts[features.selected,:])
    Y = ExpressionMatrix(convert(Matrix{Float32}, counts), featurenames)
    # make cell sketches
    T = n_bits ==  64 ? BitVec64  :
        n_bits == 128 ? BitVec128 :
        n_bits == 256 ? BitVec256 :
        n_bits == 512 ? BitVec512 :
        @assert(false)
    if transformer == :log1p
        transformer = LogT(1.0f0)
    else
        @assert(transformer == :ftt)
        transformer = FTT()
    end
    # preprocess expression profiles
    preproc = Preprocessor(Y, transformer, PCA(n_dims, randomize=randomize), normalize, standardize, Float32(scalefactor))
    X = preprocess(preproc, Y, false)
    # hash preprocessed data
    lshashes = LSHash{T}[]
    for _ in 1:n_lshashes
        P = generate_random_projections(n_bits, n_dims, superbit)
        Z = Vector{T}(undef, n)
        sketch!(Z, X, P)
        push!(lshashes, LSHash(P, HammingIndexes.HammingIndex(Z, index=index)))
    end
    # compress counts if required
    if keep_counts
        ccounts = CMatrix(counts, mmap=true)
    else
        ccounts = nothing
    end
    # create a cell index
    return CellIndex(preproc, lshashes, metadata, ccounts, RuntimeStats())
end

"""
    save(filename::AbstractString, index::CellIndex)

Save `index` to a file.

See also `load`.
"""
function save(filename::AbstractString, index::CellIndex)
    Base.open(filename, "w") do output
        serialize(output, index.preproc)
        serialize(output, index.lshashes)
        serialize(output, index.metadata)
        serialize(output, index.rtstats)
        if index.counts !== nothing
            counts = index.counts
            write(output, Int64(counts.size[1]))
            write(output, Int64(counts.size[2]))
            write(output, Int64(sizeof(counts.data)))
            write(output, counts.data)
        end
    end
    return nothing
end

"""
    load(filename::AbstractString; mmap::Bool=true)

Load an index from `filename`.

If `mmap` is `true`, the expression matrix stored in the file will be loaded
using memory-mapped file to save the memory space.

See also `save`.
"""
function load(filename::AbstractString; mmap::Bool=true)
    Base.open(filename, "r") do input
        preproc = deserialize(input)
        lshashes = deserialize(input)
        metadata = deserialize(input)
        rtstats = deserialize(input)
        if eof(input)
            counts = nothing
        else
            m = read(input, Int64)
            n = read(input, Int64)
            len = read(input, Int64)
            if mmap
                data = Mmap.mmap(filename, Vector{UInt8}, len, position(input), grow=false)
            else
                data = read(input, len)
            end
            counts = CMatrix{Int32}(data, (m, n))
        end
        return CellIndex(preproc, lshashes, metadata, counts, rtstats)
    end
end

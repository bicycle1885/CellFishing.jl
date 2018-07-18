module CellFishing

using CodecZlib: GzipDecompressorStream
using CodecZstd: ZstdDecompressorStream

if VERSION > v"0.7-"
    using Arpack: svds
    using LinearAlgebra: lu!, qr!, svd
    using Random: MersenneTwister
    using Serialization: serialize, deserialize
    using SparseArrays: SparseMatrixCSC
    using Statistics: mean, std
    macro f(ex) esc(ex) end
else
    using Compat: undef, minimum, maximum, sum, mean, std, findfirst
    const lu! = lufact!
    const qr! = qrfact!
    macro f(ex)
        function rec(ex)
            if ex isa Expr
                if ex.head == :.
                    Expr(:ref, rec(ex.args[1]), ex.args[2])
                else
                    Expr(ex.head, rec.(ex.args)...)
                end
            else
                ex
            end
        end
        esc(rec(ex))
    end
end

include("svd.jl")
include("bitvectors.jl")
include("HammingIndexes.jl")

const HammingIndex = HammingIndexes.HammingIndex
const NearestNeighbors = HammingIndexes.NearestNeighbors
const bitsof = HammingIndexes.bitsof
const prefetch = HammingIndexes.prefetch

# Compute Hamming distance between x and y.
hammdist(x::T, y::T) where {T<:BitVec} = count_ones(x ⊻ y)

# Compute approximated angle between x and y.
approxangle(x::T, y::T) where {T<:BitVec} = hammdist(x, y) * Float32(pi) / bitsof(T)


# Expression Matrix
# -----------------

# A matrix wrapper.
struct ExpressionMatrix{T}
    data::AbstractMatrix{T}
    featuremap::Dict{String,Int}
    featurenames::Vector{String}

    function ExpressionMatrix(data::AbstractMatrix{T}, featurenames) where T
        m, n = size(data)
        if m != length(featurenames)
            throw(ArgumentError("wrong number of feature names"))
        end
        featuremap = Dict{String,Int}()
        for (i, feature) in enumerate(featurenames)
            if haskey(featuremap, feature)
                throw(ArgumentError("found duplicated feature names: '$(feature)'"))
            end
            featuremap[feature] = i
        end
        return new{T}(data, featuremap, featurenames)
    end
end

Base.size(M::ExpressionMatrix) = size(M.data)
Base.size(M::ExpressionMatrix, d::Integer) = size(M.data, d)
Base.convert(::Type{ExpressionMatrix{S}}, M::ExpressionMatrix{T}) where {S,T} =
    ExpressionMatrix(convert(Matrix{S}, M.data), M.featurenames)


# Transformation
# --------------

# Logarithmic transformtion
struct LogT
    pseudocount::Float32
end
function transform!(t::LogT, X::Matrix{Float32})
    @fastmath @inbounds for i in 1:length(X)
        X[i] = log(X[i] + t.pseudocount)
    end
    return X
end

# Freeman-Tukey transformation
struct FTT end
function transform!(::FTT, X::Matrix{Float32})
    @fastmath @inbounds for i in 1:length(X)
        X[i] = sqrt(X[i]) + sqrt(X[i]+1.0f0)
    end
    return X
end


# Dimensionality Reduction
# ------------------------

mutable struct PCA
    dims::Int
    randomize::Bool
    proj::Matrix{Float32}
    PCA(dims::Integer; randomize::Bool=true) = new(dims, randomize)
end

# NOTE: This assumes features are already centered (if needed).
function fit!(pca::PCA, X::AbstractMatrix)
    if pca.randomize
        pca.proj, = rsvd(X, pca.dims)
    else
        pca.proj, = tsvd(X, pca.dims)
    end
    return pca
end
reducedims(pca::PCA, X::AbstractMatrix) = pca.proj'X


# Preprocessor
# ------------

# permutation, log transformation, scaling, projection, etc.
struct Preprocessor
    featurenames::Vector{String}
    transformer::Union{LogT,FTT}
    dimreducer::PCA
    normalize::Bool
    standardize::Bool
    scalefactor::Float32     # used for normalization
    mean::Vector{Float32}    # used for PCA
    invstd::Vector{Float32}  # used for standardization
end

function Preprocessor(
        Y::ExpressionMatrix,
        transformer,
        dimreducer,
        normalize::Bool,
        standardize::Bool,
        scalefactor::Float32)
    @assert Y.data isa Matrix{Float32}
    X::Matrix{Float32} = copy(Y.data)
    # normalize cell counts
    if normalize
        X .*= scalefactor ./ sum(X, dims=1)
    end
    # apply variance-stabilizing transformation
    transform!(transformer, X)
    # standardize features
    μ = vec(mean(X, dims=2))
    X .-= μ
    if standardize
        σ = vec(std(X, dims=2))
        if any(σ .== 0)
            throw(ArgumentError("found $(sum(σ .== 0)) features with no variance; filter out these features"))
        end
        invstd = inv.(σ)
        X .*= invstd
    else
        invstd = Float32[]
    end
    # fit dimreducer
    fit!(dimreducer, X)
    return Preprocessor(
        Y.featurenames,
        transformer,
        dimreducer,
        normalize,
        standardize,
        scalefactor,
        μ, invstd)
end

function preprocess(proc::Preprocessor, Y::ExpressionMatrix, inferparams::Bool)
    perm = zeros(Int, length(proc.featurenames))
    for (i, name) in enumerate(proc.featurenames)
        perm[i] = Y.featuremap[name]
    end
    X = permuterows(perm, Y.data)::Matrix{Float32}  # should be inferable
    if proc.normalize
        X .*= proc.scalefactor ./ sum(X, dims=1)
    end
    transform!(proc.transformer, X)
    μ = inferparams ? vec(mean(X, dims=2)) : proc.mean
    if proc.standardize
        if inferparams
            σ = vec(std(X, dims=2))
            for i in 1:length(σ)
                σ[i] = ifelse(σ[i] > 0, σ[i], proc.invstd[i])
            end
            invstd = inv.(σ)
        else
            invstd = proc.invstd
        end
        @inbounds for j in 1:size(X, 2), i in 1:size(X, 1)
            X[i,j] = (X[i,j] - μ[i]) * invstd[i]
        end
    else
        X .-= μ
    end
    return reducedims(proc.dimreducer, X)
end

indims(p::Preprocessor) = length(p.featurenames)
outdims(p::Preprocessor) = p.dimreducer.dims

permuterows(perm::Vector{Int}, Y::AbstractMatrix) = convert(Matrix{Float32}, Y[perm,:])

function permuterows(perm::Vector{Int}, Y::Matrix{Float32})
    m = length(perm)
    n = size(Y, 2)
    Y_permuted = Matrix{Float32}(undef, m, n)
    @inbounds for j in 1:n
        L = 16
        col = (j-1) * size(Y, 1)
        for i in 1:m-L
            prefetch(pointer(Y, perm[i+L] + col))
            Y_permuted[i,j] = Y[perm[i] + col]
        end
        for i in m-L+1:m
            Y_permuted[i,j] = Y[perm[i] + col]
        end
    end
    return Y_permuted
end


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
    hammingindex::HammingIndex{T}
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
        push!(Ps, Matrix(@f qr!(M).Q)')
    end
    if r > 0
        M = randn(Float32, D, r)
        push!(Ps, Matrix(@f qr!(M).Q)')
    end
    return vcat(Ps...)
end


# Feature Selection
# -----------------

"""
A list of features.

Fields
------

- `names`: feature names.
- `selected`: a bit vector of selected features.

`names` and `selected` must have the same number of elements.
"""
struct Features
    names::Vector{String}
    selected::BitVector

    function Features(names, selected)
        if length(names) != length(selected)
            throw(ArgumentError("mismatching length"))
        end
        return new(names, selected)
    end
end

"""
    Features(names)

Create a features object from feature names `names`.

All features are not selected by default. To select and unselect some features,
use `addfeatures!` and `dropfeatures!`, respectively.

`selectfeatures` algorithmically selects important features from a count matrix.
"""
Features(names::AbstractVector{String}) = Features(names, falses(length(names)))

Base.copy(features::Features) = Features(copy(features.names), copy(features.selected))

function Base.show(io::IO, features::Features)
    print(io, summary(features), "(<n-features=$(nfeatures(features)),n-selected=$(nselected(features))>)")
end

nfeatures(features::Features) = length(features.names)
nselected(features::Features) = sum(features.selected)
selectedfeatures(features::Features) = features.names[features.selected]

"""
    addfeatures(features, names) -> Features

Create a new features object by adding `names` to `features`.

See also `dropfeatures` and `addfeatures!`.
"""
addfeatures(features::Features, featurenames::AbstractVector{String}) = addfeatures!(copy(features), featurenames)

"""
    addfeatures!(features, names) -> Features

Add `names` to `features` in place.

See also `dropfeatures!` and `addfeatures`.
"""
function addfeatures!(features::Features, featurenames::AbstractVector{String})
    for name in featurenames
        i = findfirst(isequal(name), features.names)
        if i == nothing
            throw(ArgumentError("not found feature '$(name)'"))
        end
        features.selected[i] = true
    end
    return features
end

"""
    dropfeatures(features, names) -> Features

Create a new features object by dropping `names` from `features`.

See also `addfeatures` and `dropfeatures!`.
"""
dropfeatures(features::Features, featurenames::AbstractVector{String}) = dropfeatures!(copy(features), featurenames)


"""
    dropfeatures!(features, names) -> Features

Drop `names` from `features` in place.

See also `addfeatures!` and `dropfeatures`.
"""
function dropfeatures!(features::Features, featurenames::AbstractVector{String})
    for name in featurenames
        i = findfirst(isequal(name), features.names)
        if i == nothing
            throw(ArgumentError("not found feature '$(name)'"))
        end
        features.selected[i] = false
    end
    return features
end

"""
    selectfeatures(counts, featurenames; n_min_features=cld(size(counts, 1), 10)) -> Features

Select features from `counts`.

Arguments
---------

- `counts`: transcriptome expression matrix (features x cells).
- `featurenames`: vector of feature names.
- `n_min_features`: number of minimum features [default: 10% of all features].
"""
function selectfeatures(counts::AbstractMatrix,
                        featurenames::AbstractVector{String};
                        n_min_features::Real=cld(size(counts, 1), 10))
    M = size(counts, 1)
    if M != length(featurenames)
        throw(ArgumentError("mismatching featurenames size"))
    elseif !(1 ≤ n_min_features ≤ M)
        throw(ArgumentError("invalid n_min_features"))
    end
    mincounts, maxcounts = feature_extrema(counts)
    lowerbound = sort(maxcounts, rev=true)[n_min_features]
    return Features(featurenames, (maxcounts .≥ lowerbound) .& (maxcounts .> mincounts))
end

function feature_extrema(counts::AbstractMatrix)
    return vec(minimum(counts, dims=2)), vec(maximum(counts, dims=2))
end

# https://github.com/JuliaLang/julia/issues/27836
function feature_extrema(counts::SparseMatrixCSC{T}) where T <: Real
    m, n = size(counts)
    mincounts = fill(typemax(T), m)
    maxcounts = fill(typemin(T), m)
    for col in 1:n
        for j in counts.colptr[col]:counts.colptr[col+1]-1
            i = counts.rowval[j]
            val = counts.nzval[j]
            mincounts[i] = min(val, mincounts[i])
            maxcounts[i] = max(val, maxcounts[i])
        end
    end
    return mincounts, maxcounts
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
    # run-time statistics
    rtstats::RuntimeStats
end

nbits(index::CellIndex) = bitsof(bitvectype(index.lshashes[1]))
ncells(index::CellIndex) = length(index.lshashes[1].hammingindex)

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
        index::Bool=true,)
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
    Y = ExpressionMatrix(convert(Matrix{Float32}, Matrix(counts[features.selected,:])), featurenames)
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
        push!(lshashes, LSHash(P, HammingIndex(Z, index=index)))
    end
    # create a cell index
    return CellIndex(preproc, lshashes, metadata, RuntimeStats())
end

"""
    save(filename::AbstractString, index::CellIndex)

Save `index` to a file.

See also `load`.
"""
function save(filename::AbstractString, index::CellIndex)
    open(io->serialize(io, index), filename, "w")
    return nothing
end

"""
    load(filename::AbstractString)

Load an index from `filename`.

The file may be compressed with gzip or zstd; the compression format will be
detected by the extension of the file.

See also `save`.
"""
function load(filename::AbstractString)
    open(filename, "r") do file
        if endswith(filename, ".gz")
            stream = GzipDecompressorStream(file)
        elseif endswith(filename, ".zst")
            stream = ZstdDecompressorStream(file)
        else
            stream = file
        end
        return deserialize(stream)
    end
end


# k-NN Searcher
# -------------

"""
A set of k-nearest neighboring cells.

Fields
------

- `indexes`: cell indexes (k x n-queries).
- `hammingdistances`: Hamming distances (k x n-queries).
"""
struct NearestCells
    indexes::Matrix{Int}
    hammingdistances::Matrix{Int16}
    # some scores?

    NearestCells(k, n) = new(zeros(Int, k, n), zeros(Int16, k, n))
end

"""
    findneighbors(
        k::Integer,
        counts::AbstractMatrix,
        featurenames::AbstractVector{String},
        index::CellIndex;
        inferparams::Bool=true
    ) -> NearestCells

Find `k`-nearest neighboring cells from `index`.

If `inferparams=true`, feature (gene) parameters are inferred from the `counts`
argument.  Note that `counts` should not be biased and have enough cells to
properly infer the parameters.
"""
function findneighbors(
        k::Integer,
        counts::AbstractMatrix,
        featurenames::AbstractVector{String},
        index::CellIndex;
        inferparams::Bool=true)
    return findneighbors(k, ExpressionMatrix(counts, featurenames), index; inferparams=inferparams)
end

function findneighbors(k::Integer, Y::ExpressionMatrix, index::CellIndex; inferparams::Bool=true)
    if k < 0
        throw(ArgumentError("negative k"))
    end
    n = size(Y, 2)
    L = length(index.lshashes)
    T = bitvectype(first(index.lshashes))
    @assert L ≥ 1
    rtstats = index.rtstats
    tic!(rtstats)
    # preprocess
    X = preprocess(index.preproc, Y, inferparams)
    # allocate temporary memories
    neighbors = Matrix{Int}(undef, k * L, n)
    neighbors_tmp = Matrix{Int}(undef, k, n)
    Z = Matrix{T}(undef, L, n)
    Z_tmp = Vector{T}(undef, n)
    for l in 1:L
        rtstats.preproc_time += tic!(rtstats)
        findneighbors!(neighbors_tmp, Z_tmp, X, index.lshashes[l])
        neighbors[k*(l-1)+1:k*l,:] = neighbors_tmp
        Z[l,:] = Z_tmp
        rtstats.search_time += tic!(rtstats)
    end
    @assert all(neighbors .> 0)
    # Rank neighbor cells.
    ncs = rankcells!(neighbors, Z, index.lshashes)
    rtstats.rank_time += tic!(rtstats)
    return ncs
end

function rankcells!(neighbors::Matrix{Int}, Z::Matrix{T}, lshashes::Vector{LSHash{T}}) where T
    L, N = size(Z)
    k = div(size(neighbors, 1), L)
    @assert size(neighbors) == (k * L, N)
    ncs = NearestCells(k, N)
    fill!(ncs.hammingdistances, L * bitsof(T) + 1)
    for j in 1:N
        last = 0
        sort!(@view neighbors[:,j])
        for i in 1:k*L
            x = neighbors[i,j]
            if x == last
                continue
            end
            d = 0
            for l in 1:L
                d += hammdist(Z[l,j], lshashes[l].hammingindex[x])
            end
            if ncs.indexes[k,j] == 0 || d < ncs.hammingdistances[k,j]
                ncs.indexes[k,j] = x
                ncs.hammingdistances[k,j] = d
                l = k - 1
                while l ≥ 1 && ncs.hammingdistances[l,j] > ncs.hammingdistances[l+1,j]
                    ncs.indexes[l,j], ncs.indexes[l+1,j] = ncs.indexes[l+1,j], ncs.indexes[l,j]
                    ncs.hammingdistances[l,j], ncs.hammingdistances[l+1,j] = ncs.hammingdistances[l+1,j], ncs.hammingdistances[l,j]
                    l -= 1
                end
            end
            last = x
        end
        #@assert issorted(ncs.hammingdistances[:,j])
    end
    return ncs
end

end

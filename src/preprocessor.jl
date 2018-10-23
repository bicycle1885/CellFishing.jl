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

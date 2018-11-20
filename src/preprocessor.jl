# Preprocessor
# ============

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
    n::Int
    featurenames::Vector{String}
    sums::Vector{Float32}
    soss::Vector{Float32}
    transformer::Union{LogT,FTT}
    dimreducer::PCA
    normalize::Bool
    standardize::Bool
    scalefactor::Float32
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
    m, n = size(X)
    sums = vec(sum(X, dims=2))
    soss = vec(sum(X.^2, dims=2))
    μ = sums ./ n
    X .-= μ
    if standardize
        σ = sqrt.((soss .- 2 .* μ .* sums) ./ n .+ μ.^2)
        if any(σ .== 0)
            throw(ArgumentError("found $(sum(σ .== 0)) features with no variance; filter out these features"))
        end
        X .*= inv.(σ)
    end
    # fit dimreducer
    fit!(dimreducer, X)
    return Preprocessor(
        n,
        Y.featurenames,
        sums,
        soss,
        transformer,
        dimreducer,
        normalize,
        standardize,
        scalefactor,
    )
end

function preprocess(proc::Preprocessor, Y::ExpressionMatrix, inferstats::Symbol)
    perm = zeros(Int, length(proc.featurenames))
    for (i, name) in enumerate(proc.featurenames)
        perm[i] = get(Y.featuremap, name, 0)
    end
    X = permuterows(perm, Y.data)::Matrix{Float32}  # should be inferable
    if proc.normalize
        X .*= proc.scalefactor ./ sum(X, dims=1)
    end
    transform!(proc.transformer, X)
    m, n = size(X)
    @assert length(proc.sums) == length(proc.soss) == m
    if inferstats == :query
        μ = vec(mean(X, dims=2))
        σ = vec(std(X, dims=2))
        @inbounds for i in 1:m
            # insert the std of the database if no variability in the query cells
            if σ[i] == 0
                σ[i] = sqrt((proc.soss[i] - 2 * μ[i] * proc.sums[i]) / proc.n + μ[i]^2)
            end
        end
    elseif inferstats == :database
        μ = proc.sums ./ proc.n
        σ = sqrt.((proc.soss .- 2 .* μ .* proc.sums) ./ proc.n .+ μ.^2)
    else
        @assert inferstats == :both
        sums = vec(sum(X, dims=2))
        soss = vec(sum(X.^2, dims=2))
        μ = zeros(Float32, m)
        σ = zeros(Float32, m)
        @inbounds for i in 1:m
            μ[i], σ[i] = mean_and_std(
                proc.sums[i], proc.soss[i], proc.n,
                sums[i], soss[i], n,
            )
        end
    end
    if proc.standardize
        invstd = inv.(σ)
        @inbounds for j in 1:n, i in 1:m
            X[i,j] = (X[i,j] - μ[i]) * invstd[i]
        end
    else
        X .-= μ
    end
    return reducedims(proc.dimreducer, X)
end

# sos: sum of squares
@inline function mean_and_std(
        sum1::Float32, sos1::Float32, n1::Integer,
        sum2::Float32, sos2::Float32, n2::Integer,
    )
    sum = sum1 + sum2
    d = inv(n1 + n2)
    μ = sum * d
    σ = sqrt((sos1 + sos2 - 2 * μ * sum) * d + μ^2)
    return μ, σ
end

indims(p::Preprocessor) = length(p.featurenames)
outdims(p::Preprocessor) = p.dimreducer.dims

function permuterows(perm::Vector{Int}, Y::AbstractMatrix)
    m = length(perm)
    n = size(Y, 2)
    Y_permuted = Matrix{Float32}(undef, m, n)
    @inbounds for j in 1:n, i in 1:m
        Y_permuted[i,j] = perm[i] > 0 ? Y[perm[i],j] : 0
    end
    return Y_permuted
end

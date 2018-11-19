# Differentially Expressed Gene Analysis
# ======================================

struct DEGenes
    names::Vector{String}
    means::Matrix{Float64}
    negatives::Matrix{Float64}
    positives::Matrix{Float64}
end

const DEFAULT_NUM_NEIGHBORS = 10

function finddegs(
        counts::AbstractMatrix,
        featurenames::AbstractVector{String},
        cellindexes::AbstractVector{Int},
        index::CellIndex;
        k::Integer=DEFAULT_NUM_NEIGHBORS,)
    return finddegs(ExpressionMatrix(counts, featurenames), cellindexes, index, k=k)
end

function finddegs(
        Y::ExpressionMatrix,
        cellindexes::AbstractVector{Int},
        index::CellIndex;
        k::Integer=DEFAULT_NUM_NEIGHBORS,)
    n = size(Y, 2)
    if n != length(cellindexes)
        throw(ArgumentError("mismatching size"))
    elseif index.counts === nothing
        throw(ArgumentError("no counts data"))
    elseif k < 0
        throw(ArgumentError("negative k"))
    end
    neighbors = findneighbors(k, cellindexes, index)
    featurenames = index.preproc.featurenames
    m = length(featurenames)
    # an extremely weak prior (Gamma distribution)
    α₀, β₀ = 1e-1, 1e-3  # mean=100.0, var=100000.0
    #α₀, β₀ = 1.0, 1e-2  # mean=100.0, var=10000.0
    means = zeros(m, n)
    negatives = zeros(m, n)
    positives = zeros(m, n)
    for j in 1:n
        counts_j = Y[featurenames,j]
        counts_nns = index.counts[:,neighbors.indexes[:,j]]
        counts_nns_normalized = sum(counts_j) * (counts_nns ./ sum(counts_nns, dims=1))
        counts_nns_minimum = minimum(counts_nns_normalized, dims=2)
        counts_nns_maximum = maximum(counts_nns_normalized, dims=2)
        α = α₀ .+ sum(counts_nns_normalized, dims=2)
        β = β₀  + k
        for i in 1:m
            nb = NegativeBinomial(α[i], 1-inv(β+1))
            #@assert mean(nb) ≈ α[i] / β
            y = counts_j[i]
            means[i,j] = α[i] / β
            negatives[i,j] = logcdf(nb, y)
            positives[i,j] = logccdf(nb, y-1)
        end
    end
    # change the base of logarithm from e (≈2.718) to 10
    negatives .*= inv(log(10))
    positives .*= inv(log(10))
    return DEGenes(copy(featurenames), means, negatives, positives)
end

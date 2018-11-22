# Differentially Expressed Gene Analysis
# ======================================

"""
A list of probabilities of genes.

Fields
------

- `names`: gene names.
- `means`: estimated means (#genes x #queries).
- `negatives`: logarithmic probabilities (the base is 10; #genes × #queries).
- `positives`: logarithmic probabilities (the base is 10; #genes × #queries).
"""
struct DEGenes
    names::Vector{String}
    means::Matrix{Float64}
    negatives::Matrix{Float64}
    positives::Matrix{Float64}
end

const DEFAULT_NUM_NEIGHBORS = 10

"""
    finddegs(
        counts::AbstractMatrix,
        featurenames::AbstractVector{<:AbstractString},
        cellindexes::AbstractVector{<:Integer},
        index::CellIndex;
        k::Integer=$(DEFAULT_NUM_NEIGHBORS),
    ) -> DEGenes

Calculate probabilities of observing the counts or extremes.

This function returns a `DEGenes` object that contains probabilities observing
the counts stored in `counts` or their extremes. The expected expression level
for each gene is estimated from the nearby cells of `cellindexes` in the
database.  That is, the lower the probability is, the infrequent observing the
count by chance is.

A `DEGenes` object has two fields to store the probabilities: `negatives` and
`positives`.  If a value of `negatives` is low, that means that the
corresponding count is supposed to be negatively regulated. Likewise, if a
value of `positives` is low, that means that the corresponding count is
supposed to be positively regulated. Note that the probabilities are
log-transformed (the base of logarithm is 10) because they may possibly be very
close to zero. For example, a probability of 0.0001 (= 0.01%) becomes -4 after
the log transformation.

Arguments
---------

Required arguments (positional):
- `counts`:
    A transcriptome expression matrix (#features × #cells).
- `featurenames`:
    A vector of feature names (#features).
    Each of the elements must correspond to a row of the `counts` matrix.
- `cellindexes`:
    A vector of cell indexes in the database to be compared with the `counts`
    matrix (#cells).
- `index`:
    A database object.
    It must store raw counts of the database cells in it, which is activated by
    passing `keep_counts=true` when creating the database object (see
    [`CellIndex`](@ref)).

Other parameters:
- `k`: The number of nearest neighbors used to estimate the mean expression level.
"""
function finddegs(
        counts::AbstractMatrix,
        featurenames::AbstractVector{<:AbstractString},
        cellindexes::AbstractVector{<:Integer},
        index::CellIndex;
        k::Integer=DEFAULT_NUM_NEIGHBORS,)
    return finddegs(ExpressionMatrix(counts, featurenames), cellindexes, index, k=k)
end

function finddegs(
        Y::ExpressionMatrix,
        cellindexes::AbstractVector{<:Integer},
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

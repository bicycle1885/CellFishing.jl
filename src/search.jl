# k-NN Searcher
# =============

"""
A set of k-nearest neighboring cells.

Fields
------

- `indexes`: cell indexes (k x #queries).
- `hammingdistances`: Hamming distances (k x #queries).
"""
struct NearestCells
    indexes::Matrix{Int}
    hammingdistances::Matrix{Int16}
    # some scores?

    NearestCells(k, n) = new(zeros(Int, k, n), zeros(Int16, k, n))
end

"""
    similarities(ncs::NearestCells, index::CellIndex) -> Matrix{Float32}

Estimate cosine similarities between queries and neighbors from the Hamming distance.
"""
function similarities(ncs::NearestCells, index::CellIndex)
    T = nbits(index) * length(index.lshashes)
    return cos.(ncs.hammingdistances .* (Float32(pi) ./ T))
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

function findneighbors(k::Integer, I::AbstractVector{Int}, index::CellIndex)
    if k < 0
        throw(ArgumentError("negative k"))
    end
    n = length(I)
    L = length(index.lshashes)
    @assert L ≥ 1
    T = CellFishing.bitvectype(first(index.lshashes))
    neighbors = zeros(Int, k * L, n)
    Z = Matrix{T}(undef, L, n)
    nns = HammingIndexes.NearestNeighbors(k)
    for l in 1:L
        lshash = index.lshashes[l]
        for j in 1:n
            h = lshash.hammingindex[I[j]]
            Z[l,j] = h
            HammingIndexes.findknn!(nns, h, lshash.hammingindex)
            for i in 1:k
                neighbors[k*(l-1)+i,j] = nns.indexes[i]
            end
        end
    end
    @assert all(neighbors .> 0)
    return CellFishing.rankcells!(neighbors, Z, index.lshashes)
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

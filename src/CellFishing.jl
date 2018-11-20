module CellFishing

using Arpack: svds
using LinearAlgebra: lu!, qr!, svd
using Random: MersenneTwister
using Serialization: serialize, deserialize
using SparseArrays: SparseMatrixCSC
using Statistics: mean, std
using Mmap: Mmap
using Blosc: libblosc
using Distributions: NegativeBinomial, logcdf, logccdf

include("svd.jl")
include("bitvectors.jl")
include("HammingIndexes.jl")
include("cmatrix.jl")
include("ematrix.jl")
include("preprocessor.jl")
include("features.jl")
include("index.jl")
include("search.jl")
include("degenes.jl")

using .HammingIndexes: bitsof, prefetch

# Compute Hamming distance between x and y.
hammdist(x::T, y::T) where {T<:BitVec} = count_ones(x ⊻ y)

# Compute approximated angle between x and y.
approxangle(x::T, y::T) where {T<:BitVec} = hammdist(x, y) * Float32(pi) / bitsof(T)

end  # module CellFishing

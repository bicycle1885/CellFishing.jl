module CellFishing

using Arpack: svds
using LinearAlgebra: lu!, qr!, svd
using Random: MersenneTwister
using Serialization: serialize, deserialize
using SparseArrays: SparseMatrixCSC
using Statistics: mean, std
using Mmap: Mmap
using CodecZlib: GzipDecompressorStream
using CodecZstd: ZstdDecompressorStream
using Blosc: libblosc
using Distributions: NegativeBinomial, logcdf, logccdf

include("svd.jl")
include("bitvectors.jl")
include("HammingIndexes.jl")
include("ematrix.jl")
include("preprocessor.jl")
include("features.jl")
include("cmatrix.jl")
include("index.jl")
include("search.jl")
include("degenes.jl")

const bitsof = HammingIndexes.bitsof
const prefetch = HammingIndexes.prefetch

# Compute Hamming distance between x and y.
hammdist(x::T, y::T) where {T<:BitVec} = count_ones(x ⊻ y)

# Compute approximated angle between x and y.
approxangle(x::T, y::T) where {T<:BitVec} = hammdist(x, y) * Float32(pi) / bitsof(T)


end  # module CellFishing

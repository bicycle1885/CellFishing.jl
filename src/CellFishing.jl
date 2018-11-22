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

using .HammingIndexes: bitsof

end  # module CellFishing

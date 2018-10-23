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

function Base.getindex(M::ExpressionMatrix, featurenames::AbstractVector{<:AbstractString}, j::Integer)
    m = length(featurenames)
    v = Vector{eltype(M)}(undef, m)
    for i in 1:m
        v[i] = M.data[M.featuremap[featurenames[i]],j]
    end
    return v
end

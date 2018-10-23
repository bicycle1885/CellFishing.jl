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
    print(io, summary(features), "(<#features=$(nfeatures(features)),#selected=$(nselected(features))>)")
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

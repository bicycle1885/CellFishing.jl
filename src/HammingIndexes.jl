module HammingIndexes

using Compat: undef, copyto!, Nothing


# Simple List
# -----------

mutable struct List{T}
    data::Vector{T}
    size::Int

    function List{T}(sizehint::Integer) where T
        data = Vector{T}(undef, sizehint)
        return new(data, 0)
    end
end

Base.length(list::List) = list.size
Base.start(list::List) = 1
Base.done(list::List, i::Int) = i > list.size
Base.next(list::List, i::Int) = list.data[i], i + 1

function Base.empty!(list::List)
    list.size = 0
    return list
end

# This is fast due to inlining.
@inline function Base.push!(list::List{T}, val) where T
    if list.size ≥ length(list.data)
        resize!(list.data, length(list.data) * 2)
    end
    @inbounds list.data[list.size+1] = val
    list.size += 1
    return list
end


# Bit Block
# ---------

struct Block
    large::UInt32
    smalls::NTuple{4,UInt8}
    chunks::NTuple{4,UInt64}
end

function Block(chunks::NTuple{4,UInt64}, offset::Int)
    a =     convert(UInt8, count_ones(chunks[1]))
    b = a + convert(UInt8, count_ones(chunks[2]))
    c = b + convert(UInt8, count_ones(chunks[3]))
    Block(offset & ~UInt32(0), (UInt8(offset >>> 32), a, b, c), chunks)
end

const BITS_PER_CHUNK =  64
const BITS_PER_BLOCK = 256


# Succinct Bit Vector
# -------------------

struct SucVector <: AbstractVector{Bool}
    blocks::Vector{Block}
    len::Int
end

function SucVector(bv::BitVector)
    len = length(bv)
    @assert len ≤ 2^40
    n_blocks = cld(len, BITS_PER_BLOCK)
    blocks = Vector{Block}(undef, n_blocks)
    offset = 0
    for i in 1:n_blocks
        chunks = read_4chunks(bv, (i - 1) * BITS_PER_BLOCK + 1)
        blocks[i] = Block(chunks, offset)
        for j in 1:4
            offset += count_ones(chunks[j])
        end
    end
    return SucVector(blocks, len)
end

function read_4chunks(src::BitVector, from::Int)
    @assert rem(from, BITS_PER_BLOCK) == 1
    # NOTE: this method depends on the internal data layout of BitVector
    # (but much faster).
    if length(src) >= from + BITS_PER_BLOCK - 1
        i = div(from - 1, BITS_PER_CHUNK) + 1
        a = src.chunks[i]
        b = src.chunks[i+1]
        c = src.chunks[i+2]
        d = src.chunks[i+3]
    else
        a = read_chunk(src, from)
        b = read_chunk(src, from + BITS_PER_CHUNK * 1)
        c = read_chunk(src, from + BITS_PER_CHUNK * 2)
        d = read_chunk(src, from + BITS_PER_CHUNK * 3)
    end
    return a, b, c, d
end

function read_chunk(src, from::Int)
    @assert BITS_PER_CHUNK == sizeof(UInt64) * 8
    chunk = UInt64(0)
    for k in 0:BITS_PER_CHUNK-1
        i = from + k
        chunk >>= 1
        if i ≤ endof(src) && src[i]
            chunk |= UInt64(1) << 63
        end
    end
    return chunk
end

Base.size(v::SucVector) = (v.len,)

@inline function block_id(i)
    j = Int(i - 1)
    (j >>> 8) + 1, (j & 0b11111111) + 1
end

@inline function chunk_id(i)
    j = Int(i - 1)
    (j >>> 6) + 1, (j & 0b00111111) + 1
end

@inline function Base.getindex(v::SucVector, i::Int)
    q, r = block_id(i)
    @inbounds block = v.blocks[q]
    q, r = chunk_id(r)
    @inbounds chunk = block.chunks[q]
    return (chunk >> (r - 1)) & 1 == 1
end

# O(1) bit counting.
@inline function rank1(bv::SucVector, i::Int)
    i = clamp(i, 0, bv.len)
    if i == 0
        return 0
    end
    q, r = block_id(i)
    @inbounds begin
        block = bv.blocks[q]
        # large block
        ret = Int(block.large) + Int(block.smalls[1]) << 32
        # small block
        q, r = chunk_id(r)
        ret += ifelse(q == 1, 0x00, block.smalls[q])
        # remaining bits
        chunk = block.chunks[q]
        ret += count_ones(chunk & ((UInt64(1) << r) - 1))
        return ret
    end
end


# SubIndex
# --------

struct SubIndex
    # substring length
    s::Int
    # direct addressing buckets
    filled::SucVector
    offsets::Vector{Int32}  # offsets of buckets
    buckets::Vector{Int32}  # concatenated buckets

    function SubIndex(s::Int, filled::SucVector, buckets::Vector{Vector{Int32}})
        @assert s ≤ 32
        @assert length(buckets) == rank1(filled, length(filled))
        offsets = Vector{Int32}(undef, length(buckets)+1)
        buckets_concat = zeros(Int32, sum(length(b) for b in buckets) + 16)  # margin for prefetching
        offset = 1
        for (i, bucket) in enumerate(buckets)
            offsets[i] = offset
            copyto!(buckets_concat, offset, bucket, 1, length(bucket))
            offset += length(bucket)
        end
        offsets[length(buckets)+1] = offset
        return new(s, filled, offsets, buckets_concat)
    end
end

bitsof(::Type{T}) where {T} = sizeof(T) * 8

function findat!(results::List{UnitRange{Int32}}, r::Int, g::UInt32, index::SubIndex)
    @assert r ≤ index.s ≤ 32
    x_max = (UInt32(1) << index.s) - UInt32(1)
    # bit vector to flip r bits of g
    x = (UInt32(1) << r) - UInt32(1)
    while true
        @assert count_ones(x) == r
        k = Int(g ⊻ x) + 1
        if index.filled[k]
            i = rank1(index.filled, k)
            @inbounds push!(results, index.offsets[i]:index.offsets[i+1]-1)
        end
        x_next = snoob(x)
        if x_next == x || x_next > x_max
            break
        end
        x = x_next
    end
    return results
end

function findwithin!(results::List{UnitRange{Int32}}, r::Int, g::UInt32, index::SubIndex)
    @assert r ≤ index.s ≤ 32
    for r′ in 0:r
        findat!(results, r′, g, index)
    end
    return results
end

# https://www.slideshare.net/gkumar007/bits-next-higher-presentation
function snoob(x::UInt32)
    y = x & -x
    if y == 0
        return x
    end
    z = x + y
    w = div(x ⊻ z, y) >> 2
    return max(z | w, x)
end


# HammingIndex
# ------------

struct HammingIndex{T}
    # data points in the Hamming space
    data::Vector{T}
    # offsets of the j-th substring
    offsets::Vector{Int}
    # indexes of substrings
    subindexes::Vector{SubIndex}
end

Base.length(index::HammingIndex) = length(index.data)
Base.getindex(index::HammingIndex, i::Integer) = index.data[i]

isindexed(index::HammingIndex) = !isempty(index.subindexes)

function HammingIndex(data::Vector{T}; index::Bool=true) where T
    if !index
        # no index
        return HammingIndex{T}(data, Int[], SubIndex[])
    end
    @assert length(data) ≤ typemax(Int32)
    @assert sizeof(T) ≥ 2
    # q-bit Hamming space
    q = bitsof(T)
    s = clamp(round(Int, log2(length(data))), 2, 32)
    m, r = divrem(q, s)
    offsets = [0]
    for j in 1:m
        if r > 0
            s_j = s + 1
            r -= 1
        else
            s_j = s
        end
        push!(offsets, offsets[end] + s_j)
    end
    subindexes = SubIndex[]
    for j in 1:m
        s_j = offsets[j+1] - offsets[j]
        offset = offsets[j]
        mask = substrmask(offsets, j)
        # first pass: determine buckets to fill
        filled = falses(2^s_j)
        for h in data
            h_j = ((h >> offset) % UInt32) & mask
            filled[Int(h_j) + 1] = true
        end
        # second pass: fill buckets
        filled_sv = SucVector(filled)
        buckets = [Int32[] for _ in 1:rank1(filled_sv, length(filled_sv))]
        for (i, h) in enumerate(data)
            h_j = ((h >> offset) % UInt32) & mask
            bucket = buckets[rank1(filled_sv, Int(h_j) + 1)]
            push!(bucket, i)
        end
        push!(subindexes, SubIndex(s_j, filled_sv, buckets))
    end
    return HammingIndex{T}(data, offsets, subindexes)
end

substrlen(offsets, j) = offsets[j+1] - offsets[j]
substrmask(offsets, j) = (UInt32(1) << substrlen(offsets, j)) - UInt32(1)


# Nearest Neighbors
# -----------------

# the result of k-NN search
mutable struct NearestNeighbors
    indexes::Vector{Int}
    distances::Vector{Int}
    n_compared::Int     # number of comparisons
    time_elapsed::UInt  # elapsed time in nanosecond
end

NearestNeighbors(k::Integer) = NearestNeighbors(zeros(Int, k), zeros(Int, k), 0, 0)
Base.length(nns::NearestNeighbors) = length(nns.indexes)

# Tidy up `nns` when the last element is updated.
function tidyup!(nns::NearestNeighbors)
    j = length(nns) - 1
    while j ≥ 1 && nns.distances[j] > nns.distances[j+1]
        nns.indexes[j], nns.indexes[j+1] = nns.indexes[j+1], nns.indexes[j]
        nns.distances[j], nns.distances[j+1] = nns.distances[j+1], nns.distances[j]
        j -= 1
    end
    return nns
end


# High-level Search Functions
# ---------------------------

function findknn(k::Int, q::T, index::HammingIndex{T}) where {T}
    if isindexed(index)
        return findknn(MultiIndexSearch(), q, index)
    else
        return findknn(LinearSearch(), q, index)
    end
end

function findknn!(nns::NearestNeighbors, q::T, index::HammingIndex{T}) where {T}
    if isindexed(index)
        return findknn!(MultiIndexSearch(), nns, q, index)
    else
        return findknn!(LinearSearch(), nns, q, index)
    end
end


# Linear Search
# -------------

struct LinearSearch end

function findwithin(::LinearSearch, r::Int, q::T, index::HammingIndex{T}) where {T}
    results = Int[]
    @inbounds for i in 1:length(index.data)
        if count_ones(q ⊻ index.data[i]) ≤ r
            push!(results, i)
        end
    end
    return results
end

function findknn(algo::LinearSearch, k::Int, q::T, index::HammingIndex{T}) where {T}
    return findknn!(algo, NearestNeighbors(k), q, index)
end

function findknn!(::LinearSearch, nns::NearestNeighbors, q::T, index::HammingIndex{T}) where {T}
    time_start = time_ns()
    @assert length(nns.indexes) == length(nns.distances)
    # initialize nearest neighbors
    fill!(nns.indexes, 0)
    fill!(nns.distances, bitsof(T) + 1)
    nns.n_compared = nns.time_elapsed = 0
    if length(nns) == 0
        @goto finish
    end
    @inbounds for i in 1:length(index.data)
        d = count_ones(q ⊻ index.data[i])
        if d < nns.distances[end]
            nns.indexes[end] = i
            nns.distances[end] = d
            tidyup!(nns)
        end
    end
    @label finish
    nns.n_compared = length(index.data)
    nns.time_elapsed = time_ns() - time_start
    return nns
end


# Multi-Index Search
# ------------------

struct MultiIndexSearch end

function findwithin(::MultiIndexSearch, r::Int, q::T, index::HammingIndex{T}) where {T}
    if r < 0
        return Int[]
    end
    # gather candidates
    m = length(index.subindexes)
    r′, a = divrem(r, m)
    results = Int[]
    ranges = List{UnitRange{Int32}}(128)
    for j in 1:m
        q_j = ((q >> index.offsets[j]) % UInt32) & substrmask(index.offsets, j)
        if j == a + 2
            r′ -= 1
        end
        subindex = index.subindexes[j]
        empty!(ranges)
        findwithin!(ranges, r′, q_j, subindex)
        for range in ranges, l in range
            push!(results, subindex.buckets[l])
        end
    end
    # filter results
    n = 0
    prev = 0
    sort!(results)
    for k in 1:length(results)
        i = results[k]
        if i != prev
            if count_ones(q ⊻ index.data[i]) ≤ r
                n += 1
                results[n] = prev = i
            end
        end
    end
    resize!(results, n)
    return results
end

function findknn(algo::MultiIndexSearch, k::Int, q::T, index::HammingIndex{T}) where {T}
    return findknn!(algo, NearestNeighbors(k), q, index)
end

function findknn!(::MultiIndexSearch, nns::NearestNeighbors, q::T, index::HammingIndex{T}) where {T}
    n_compared = 0
    time_start = time_ns()
    @assert length(nns.indexes) == length(nns.distances)
    # initialize nearest neighbors
    fill!(nns.indexes, 0)
    fill!(nns.distances, bitsof(T) + 1)
    nns.n_compared = nns.time_elapsed = 0
    if length(nns) == 0
        @goto finish
    end
    # progressively search the Hamming space for k nearest neighbors
    ranges = List{UnitRange{Int32}}(128)
    m = length(index.subindexes)
    r = r′ = a = 0
    d_max = nns.distances[end]
    while r < d_max
        q_ap1 = ((q >> index.offsets[a+1]) % UInt32) & substrmask(index.offsets, a+1)
        empty!(ranges)
        subindex = index.subindexes[a+1]
        findat!(ranges, r′, q_ap1, subindex)
        for range in ranges
            data = index.data
            buckets = subindex.buckets
            start = first(range)
            @inbounds for l in 0:length(range)-1
                prefetch(pointer(data, buckets[start+l+8]))
                i = buckets[start+l]
                p = data[i]
                d = count_ones(q ⊻ p)
                n_compared += 1
                if d < d_max
                    for j in searchsorted(nns.distances, d)
                        if nns.indexes[j] == i  # duplication
                            @goto skip
                        end
                    end
                    nns.indexes[end] = i
                    nns.distances[end] = d
                    tidyup!(nns)
                    d_max = nns.distances[end]
                    if r ≥ d_max
                        @goto finish
                    end
                end
                @label skip
            end
        end
        a += 1
        if a ≥ m
            a = 0
            r′ += 1
        end
        r += 1
        @assert r == m * r′ + a
    end
    @label finish
    nns.n_compared = n_compared
    nns.time_elapsed = time_ns() - time_start
    return nns
end

@inline function prefetch(addr::Ptr, write::Bool=false, locality::Int=0)
    ccall(
        "llvm.prefetch", llvmcall,
        Nothing, (Ptr{UInt8}, Int32, Int32, Int32),
        # address to be prefetched
        addr,
        # read mode (0), write mode (1)
        write,
        # no locality (0) - high locality (3)
        locality,
        # instruction cache (0), data cache (1)
        1,
    )
end

end

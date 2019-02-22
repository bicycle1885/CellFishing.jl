# Bit Vectors
# ===========

"""
A bit vector type to represent hash values.
"""
abstract type BitVec end

# Concrete types.
primitive type BitVec64  <: BitVec  64 end
primitive type BitVec128 <: BitVec 128 end
primitive type BitVec256 <: BitVec 256 end
primitive type BitVec512 <: BitVec 512 end

using Core: Intrinsics

const BaseInts = Union{Bool,Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64}

# Define bit operations.
for T in [:BitVec64, :BitVec128, :BitVec256, :BitVec512]
    @eval begin
        $T(x::BaseInts) = convert($T, x)
        Base.convert(::Type{$T}, x::BaseInts) = Intrinsics.zext_int($T, x)
        Base.zero(::Type{$T}) = convert($T, 0)
        Base.:&(bv1::$T, bv2::$T) = Intrinsics.and_int(bv1, bv2)
        Base.:|(bv1::$T, bv2::$T) = Intrinsics.or_int(bv1, bv2)
        Base.:⊻(bv1::$T, bv2::$T) = Intrinsics.xor_int(bv1, bv2)
        Base.:<<(bv::$T, n::BaseInts) = Intrinsics.shl_int(bv, n)
        Base.:>>(bv::$T, n::BaseInts) = Intrinsics.lshr_int(bv, n)
        Base.getindex(bv::$T, i::Integer) = Intrinsics.trunc_int(Bool, bv >> (i - 1))
        Base.rem(bv::$T, ::Type{U}) where {U<:BaseInts} = Intrinsics.trunc_int(U, bv)
        Base.count_ones(bv::$T) = Intrinsics.ctpop_int(bv) % Int
    end
end

# Intrinsics have became more stringent from Julia 1.1.
# https://github.com/JuliaLang/julia/pull/29422
Base.convert(::Type{BitVec64}, x::Union{Int64,UInt64}) = Intrinsics.bitcast(BitVec64, x)
Base.rem(bv::BitVec64, ::Type{U}) where U <: Union{Int64,UInt64} = Intrinsics.bitcast(U, bv)

function Base.rand(mt::MersenneTwister, ::Type{T}) where {T<:BitVec}
    bv::T = 0
    for _ in 1:div(sizeof(T), 8)
        bv = (bv << 64) | T(rand(mt, UInt64))
    end
    return bv
end

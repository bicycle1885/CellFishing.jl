abstract type BitVec end

primitive type BitVec32   <: BitVec   32 end
primitive type BitVec64   <: BitVec   64 end
primitive type BitVec128  <: BitVec  128 end
primitive type BitVec256  <: BitVec  256 end
primitive type BitVec512  <: BitVec  512 end
primitive type BitVec1024 <: BitVec 1024 end

const BaseInts = Union{Bool,Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64,Int128,UInt128}

using Core.Intrinsics: bitcast, zext_int, ctpop_int, trunc_int, and_int, or_int, xor_int, shl_int, lshr_int

Base.convert(::Type{BitVec32 } , x::BaseInts) = zext_int(BitVec32  , UInt32(x))
Base.convert(::Type{BitVec64 } , x::BaseInts) = zext_int(BitVec64  , UInt64(x))
Base.convert(::Type{BitVec128} , x::BaseInts) = zext_int(BitVec128 ,        x )
Base.convert(::Type{BitVec256} , x::BaseInts) = zext_int(BitVec256 ,        x )
Base.convert(::Type{BitVec512} , x::BaseInts) = zext_int(BitVec512 ,        x )
Base.convert(::Type{BitVec1024}, x::BaseInts) = zext_int(BitVec1024,        x )

BitVec32(  x::Union{Bool,BaseInts}) = convert(BitVec32  , x)
BitVec64(  x::Union{Bool,BaseInts}) = convert(BitVec64  , x)
BitVec128( x::Union{Bool,BaseInts}) = convert(BitVec128 , x)
BitVec256( x::Union{Bool,BaseInts}) = convert(BitVec256 , x)
BitVec512( x::Union{Bool,BaseInts}) = convert(BitVec512 , x)
BitVec1024(x::Union{Bool,BaseInts}) = convert(BitVec1024, x)

if VERSION > v"0.7-"
    using Random: MersenneTwister
end

Base.rand(mt::MersenneTwister, ::Type{BitVec32}) = bitcast(BitVec32, rand(mt, UInt32))
function Base.rand(mt::MersenneTwister, ::Type{T}) where {T<:BitVec}
    bv::T = 0
    for _ in 1:div(sizeof(T), 8)
        bv = (bv << 64) | T(rand(mt, UInt64))
    end
    return bv
end

for T in [BitVec32, BitVec64, BitVec128, BitVec256, BitVec512, BitVec1024]
    Base.:&(bv1::T, bv2::T) = and_int(bv1, bv2)
    Base.:|(bv1::T, bv2::T) = or_int(bv1, bv2)
    Base.:⊻(bv1::T, bv2::T) = xor_int(bv1, bv2)
    Base.:<<(bv::T, n::BaseInts) = shl_int(bv, n)
    Base.:>>(bv::T, n::BaseInts) = lshr_int(bv, n)
    Base.getindex(bv::T, i::Integer) = trunc_int(Bool, bv >> (i - 1))
    Base.rem(bv::T, ::Type{UInt32}) = trunc_int(UInt32, bv)
    if T == BitVec32
        Base.rem(bv::T, ::Type{UInt64}) = UInt64(trunc_int(UInt32, bv))
        Base.count_ones(bv::T) = Int(trunc_int(Int32, ctpop_int(bv)))
    else
        Base.rem(bv::T, ::Type{UInt64}) = trunc_int(UInt64, bv)
        Base.count_ones(bv::T) = trunc_int(Int, ctpop_int(bv))
    end
end

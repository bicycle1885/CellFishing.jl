include("../src/CellFishing.jl")

if VERSION > v"0.7-"
    using Test
    using Random
else
    using Base.Test
end

const BitVec64   = CellFishing.BitVec64
const BitVec128  = CellFishing.BitVec128
const BitVec256  = CellFishing.BitVec256
const BitVec512  = CellFishing.BitVec512
const BitVec1024 = CellFishing.BitVec1024
const HammingIndexes = CellFishing.HammingIndexes

@testset "BitVectors" begin
    for T in [BitVec64, BitVec128, BitVec256, BitVec512, BitVec1024]
        @test T <: CellFishing.BitVec
        @test convert(T, 0) isa T
        @test T(0) isa T
        @test T(false) === T(0)
        @test T(true) === T(1)
        @test (T(0b0011) & T(0b0101)) === T(0b0001)
        @test (T(0b0011) | T(0b0101)) === T(0b0111)
        @test (T(0b0011) ⊻ T(0b0101)) === T(0b0110)
        @test T(0b0010) << 1 === T(0b0100)
        @test T(0b0010) >> 1 === T(0b0001)
        @test T(0b0101)[1] === true
        @test T(0b0101)[2] === false
        @test T(0b0101) % UInt32 === UInt32(0b0101)
        @test T(0b0101) % UInt64 === UInt64(0b0101)
        @test count_ones(T(0b0101)) === 2
        @test rand(T) isa T
    end
end

@testset "HammingIndexes" begin
    linear = HammingIndexes.LinearSearch()
    mindex = HammingIndexes.MultiIndexSearch()
    srand(1234)
    data = rand(UInt64, 1000)
    db = HammingIndexes.HammingIndex(data)
    for r in 0:16
        ok = true
        for i in 1:length(data)
            x = HammingIndexes.findwithin(linear, r, data[i], db)
            y = HammingIndexes.findwithin(mindex, r, data[i], db)
            ok &= x == y
        end
        @test ok
    end
    for k in 0:16
        ok = true
        for i in 1:length(data)
            x = HammingIndexes.findknn(linear, k, data[i], db)
            y = HammingIndexes.findknn(mindex, k, data[i], db)
            # indexes may be different
            ok &= x.distances == y.distances
        end
        @test ok
    end

    srand(1234)
    data = rand(UInt64, 1_000_000)
    db = HammingIndexes.HammingIndex(data)
    queries = data[rand(1:length(data), 100)]
    for r in [16, 18, 20]
        ok = true
        for q in queries
            x = HammingIndexes.findwithin(linear, r, q, db)
            y = HammingIndexes.findwithin(mindex, r, q, db)
            ok &= x == y
        end
        @test ok
    end
    for k in 0:16
        ok = true
        for q in queries
            x = HammingIndexes.findknn(linear, k, q, db)
            y = HammingIndexes.findknn(mindex, k, q, db)
            # indexes may be different
            ok &= x.distances == y.distances
        end
        @test ok
    end
end

@testset "CellIndex" begin
    srand(1234)
    m = 100
    Y = convert(Matrix{Int32}, rand(0:1000, m, 200))
    tagnames = [string("tag:", i) for i in 1:m]
    index = CellFishing.CellIndex(
        Y,
        tagnames=tagnames,
        n_min_features=length(tagnames),
        metadata=[string("cell:", j) for j in 1:size(Y, 2)])
    @test CellFishing.nbits(index) == 128
    @test CellFishing.ncells(index) == 200
    @test index.tagnames[1] == "tag:1"
    @test index.metadata[1] == "cell:1"
    @test_throws ArgumentError CellFishing.CellIndex(Y, tagnames=["tag:1"])
    @test contains(sprint(show, index), "CellIndex")
    mktempdir() do tmpdir
        tmpfile = joinpath(tmpdir, "index")
        CellFishing.save(tmpfile, index)
        @test CellFishing.load(tmpfile) isa CellFishing.CellIndex
        run(`gzip -kq $(tmpfile)`)
        @test CellFishing.load(string(tmpfile, ".gz")) isa CellFishing.CellIndex
        run(`zstd -kq $(tmpfile)`)
        @test CellFishing.load(string(tmpfile, ".zst")) isa CellFishing.CellIndex
    end

    srand(1234)
    m = 100
    Y = convert(Matrix{Int32}, rand(0:1000, m, 200))
    tagnames = [string("tag:", i) for i in 1:m]
    n = 0
    ok = 0
    for n_bits in [64, 128, 256, 512, 1024]
        for n_dims in [64, m],
            superbit in [1, 8, 17],
            randomize in false:true,
            normalize in false:true,
            standardize in false:true,
            index in false:true
            idx = CellFishing.CellIndex(
                Y,
                tagnames=tagnames,
                n_bits=n_bits,
                n_min_features=length(tagnames),
                n_dims=n_dims, superbit=superbit, randomize=randomize,
                normalize=normalize, standardize=standardize,
                index=index)
            perm = shuffle(1:m)
            U = CellFishing.findknn(1, CellFishing.ExpressionMatrix(Y[perm,1:10], tagnames[perm]), idx)
            for i in 1:10
                j = U.indexes[1,i]
                dist = U.hammingdistances[1,i]
                @test dist == 0
                n += 1
                ok += i == j
            end
        end
    end
    println("$(ok) / $(n) = $(ok/n)")
end

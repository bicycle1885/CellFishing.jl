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
const HammingIndexes = CellFishing.HammingIndexes

@testset "BitVectors" begin
    for T in [BitVec64, BitVec128, BitVec256, BitVec512]
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

@testset "Features" begin
    m = 100
    names = string.("feature:", 1:m)
    features = CellFishing.Features(names)
    @test CellFishing.nfeatures(features) == m
    @test CellFishing.nselected(features) == 0

    # add features
    features = CellFishing.addfeatures(features, ["feature:10", "feature:32", "feature:99"])
    @test CellFishing.nfeatures(features) == m
    @test CellFishing.nselected(features) == 3
    @test CellFishing.selectedfeatures(features) == ["feature:10", "feature:32", "feature:99"]
    @test_throws ArgumentError CellFishing.addfeatures(features, ["foo", "bar"])

    # drop features
    features = CellFishing.dropfeatures(features, ["feature:32"])
    @test CellFishing.nfeatures(features) == m
    @test CellFishing.nselected(features) == 2
    @test CellFishing.selectedfeatures(features) == ["feature:10", "feature:99"]
    @test_throws ArgumentError CellFishing.dropfeatures(features, ["foo", "bar"])
end

@testset "CellIndex" begin
    srand(1234)
    m, n = 100, 200
    counts = rand(0:1000, m, n)
    featurenames = string.("feature:", 1:m)
    features = CellFishing.selectfeatures(counts, featurenames, n_min_features=m)
    index = CellFishing.CellIndex(counts, features, metadata=string.("cell:", 1:n))
    @test CellFishing.nbits(index) == 128
    @test CellFishing.ncells(index) == n
    @test index.featurenames[1] == "feature:1"
    @test index.metadata[1] == "cell:1"
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
    m, n = 100, 200
    counts = rand(0:1000, m, n)
    featurenames = string.("feature:", 1:m)
    all = ok = 0
    for n_bits in [64, 128, 256, 512],
        superbit in [1, 8, 17],
        n_dims in [64, m],
        randomize in false:true,
        normalize in false:true,
        standardize in false:true,
        index in false:true
        features = CellFishing.selectfeatures(counts, featurenames, n_min_features=m)
        idx = CellFishing.CellIndex(
            counts, features,
            n_bits=n_bits, superbit=superbit,
            n_dims=n_dims, randomize=randomize,
            normalize=normalize, standardize=standardize,
            index=index)
        perm = shuffle(1:m)
        U = CellFishing.findneighbors(1, CellFishing.ExpressionMatrix(counts[perm,1:10], featurenames[perm]), idx)
        for i in 1:10
            j = U.indexes[1,i]
            dist = U.hammingdistances[1,i]
            @test dist == 0
            all += 1
            ok += i == j
        end
    end
    @test ok / all ≈ 1.0
end

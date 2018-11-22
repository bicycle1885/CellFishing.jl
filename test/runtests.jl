using CellFishing: CellFishing, HammingIndexes
using Test
using Random: seed!, shuffle
using SparseArrays: SparseMatrixCSC

@testset "BitVectors" begin
    for T in [CellFishing.BitVec64,
              CellFishing.BitVec128,
              CellFishing.BitVec256,
              CellFishing.BitVec512]
        @test T <: CellFishing.BitVec
        @test convert(T, 0) isa T
        @test T(0) isa T
        @test zero(T) === T(0)
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
    seed!(1234)
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

    seed!(1234)
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

@testset "CMatrix" begin
    m, n = 100, 200
    seed!(1234)
    counts = rand(0:100, m, n)
    # with the default parameters
    M = CellFishing.CMatrix(counts)
    @test M isa AbstractMatrix{Int}
    @test all(M .== counts)
    for j in 1:n
        @test M[:,j] == counts[:,j]
    end
    # with specific parameters
    M = CellFishing.CMatrix(
        counts,
        level=9,
        shuffle=false,
        compressor="lz4hc",
        blocksize=m,
        nthreads=2)
    @test M isa AbstractMatrix{Int}
    @test all(M .== counts)
    for j in 1:n
        @test M[:,j] == counts[:,j]
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
    seed!(1234)
    m, n = 100, 200
    counts = rand(0:1000, m, n)
    featurenames = string.("feature:", 1:m)
    features = CellFishing.selectfeatures(counts, featurenames, n_min_features=m)
    index = CellFishing.CellIndex(counts, features, metadata=string.("cell:", 1:n))
    @test CellFishing.nbits(index) == 128
    @test CellFishing.ncells(index) == n
    @test occursin("CellFishing.CellIndex(<#cells=200, hash=128×4>)", sprint(show, index))
    #@test index.featurenames[1] == "feature:1"
    @test index.metadata[1] == "cell:1"
    @test index.counts === nothing
    @test occursin("CellIndex", sprint(show, index))
    mktempdir() do tmpdir
        tmpfile = joinpath(tmpdir, "index")
        CellFishing.save(tmpfile, index)
        @test CellFishing.load(tmpfile) isa CellFishing.CellIndex
        @test CellFishing.load(tmpfile).counts === nothing
    end
    # keep raw counts
    index = CellFishing.CellIndex(counts, features, keep_counts=true)
    @test index.counts == counts[features.selected,:]
    mktempdir() do tmpdir
        tmpfile = joinpath(tmpdir, "index")
        CellFishing.save(tmpfile, index)
        @test CellFishing.load(tmpfile) isa CellFishing.CellIndex
        @test CellFishing.load(tmpfile).counts !== nothing
    end

    seed!(12345)
    m, n = 100, 200
    counts = rand(0:1000, m, n)
    featurenames = string.("feature:", 1:m)
    n_all = n_ok = 0
    for MatType in [Matrix, SparseMatrixCSC]
        for n_bits in [64, 128, 256, 512],
            superbit in [1, 8, 17],
            n_dims in [64, m-1],
            randomize in false:true,
            normalize in false:true,
            standardize in false:true,
            index in false:true
            features = CellFishing.selectfeatures(MatType(counts), featurenames, n_min_features=m)
            idx = CellFishing.CellIndex(
                MatType(counts), features,
                n_bits=n_bits, superbit=superbit,
                n_dims=n_dims, randomize=randomize,
                normalize=normalize, standardize=standardize,
                index=index)
            perm = shuffle(1:m)
            neighbors = CellFishing.findneighbors(1, CellFishing.ExpressionMatrix(MatType(counts)[perm,:], featurenames[perm]), idx)
            totaldist::Int = 0
            for i in 1:10
                j = neighbors.indexes[1,i]
                totaldist += neighbors.hammingdistances[1,i]
                n_all += 1
                n_ok += i == j
            end
            @test totaldist ≤ 1
            S = CellFishing.similarities(neighbors, idx)
            @test all(-1 .≤ S .≤ 1)
            neighbors = CellFishing.findneighbors(3, [1, 2, 3, 5, 7], idx)
            @test neighbors.indexes[1,:] == [1, 2, 3, 5, 7]
            @test all(neighbors.hammingdistances[1,:] .== 0)
        end
    end
    @test n_ok / n_all ≈ 1.0
end

using CellFishing
using HDF5
using DocOpt

if VERSION > v"0.7-"
    using Dates
    using Random
else
    using Compat
    using Compat: @info
    struct Data
        counts::AbstractMatrix
        featurenames::Vector{String}
        cellnames::Vector{String}
    end
end

args = docopt("""
CellFishing - Cell Search Engine for scRNA-seq.

Usage:
    cellfishing build <reference>
    cellfishing search <database> <query>

Options:
    -h --help     Show this help.
    --version     Show version.
    --n-bits      The number of bits for hashing [default: 128].
    --n-lshashes  The number of locality-sensitive hashes [default: 4].
""")

function loadloom(datafile)
    HDF5.h5open(datafile) do file
        counts = Matrix(read(file, "/matrix")')
        featurenames = read(file, "/row_attrs/Gene")
        cellnames = read(file, "/col_attrs/CellID")
        @static if VERSION > v"0.7-"
            return (counts=counts, featurenames=featurenames, cellnames=cellnames)
        else
            return Data(counts, featurenames, cellnames)
        end
    end
end

function loaddata(datafile)
    if endswith(datafile, ".loom")
        return loadloom(datafile)
    else
        error("unsupported file format")
    end
end

macro message(msg, block)
    quote
        print(stderr, "  ", rpad(string($(esc(msg)), " "), 25, '\u2015'), " ")
        starttime = now()
        $(esc(block))
        println(stderr, Dates.canonicalize(Dates.CompoundPeriod(now() - starttime)))
    end
end

function writeknn(file, cellnames_query, cellnames_reference, neighbors)
    k, n = size(neighbors)
    for j in 1:n
        print(file, cellnames_query[j])
        for i in 1:k
            print(file, '\t', cellnames_reference[neighbors[i,j]])
        end
        println(file)
    end
end

srand(1234)
if args["build"]
    filename = args["<reference>"]
    println(stderr, "Build a search database from $(filename).")
    @message "Loading data" begin
        data = loaddata(filename)
    end
    @message "Selecting features" begin
        features = CellFishing.selectfeatures(data.counts, data.featurenames)
    end
    @message "Creating a database" begin
        database = CellFishing.CellIndex(data.counts, features, metadata=data.cellnames)
    end
    dbfile = string(basename(filename), ".cf")
    @message "Writing the database" begin
        CellFishing.save(dbfile, database)
    end
    println(stderr, "The serialized database is in $(dbfile).")
elseif args["search"]
    k = 10
    database = args["<database>"]
    println(stderr, "Search $(database) for $(k) neighbors.")
    @message "Loading the database" begin
        database = CellFishing.load(database)
    end
    @message "Loading query data" begin
        data = loaddata(args["<query>"])
    end
    @message "Searching the database" begin
        neighbors = CellFishing.findneighbors(k, data.counts, data.featurenames, database)
    end
    @message "Writing neighbors" begin
        writeknn(stdout, data.cellnames, database.metadata, neighbors.indexes)
    end
else
    @assert false
end

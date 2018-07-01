# CellFishing.jl ꙮ 🎣

CellFishing.jl (**cell** **fi**nder via ha**shing**) is a tool to find similar
cells of query cells based on their transcriptome expression profiles.

```julia
# Import packages.
using CellFishing
using CSV

# Load expression profiles of database cells.
data = CSV.read("database.txt", delim='\t')
cellnames = string.(names(data))
featurenames = string.(data[:,1])
counts = Matrix(data[:,2:end])

# Select features and create an index (or a database).
features = CellFishing.selectfeatures(counts, featurenames)
database = CellFishing.CellIndex(counts, features, metadata=cellnames)

# Save/load the database to/from a file (optional).
# CellFishing.save("database.cf", database)
# run(`gzip database.cf`)  # or run(`zstd database.cf`)
# database = CellFishing.load("database.cf.gz")

# Load expression profiles of query cells.
data = CSV.read("query.txt", delim='\t')
cellnames = string.(names(data))
featurenames = string.(data[:,1])
counts = Matrix(data[:,2:end])

# Search the database for similar cells; k cells will be returned per query.
k = 10
neighbors = CellFishing.findneighbors(k, counts, featurenames, database)

# Write the neighboring cells to a file.
open("neighbors.tsv", "w") do file
    println(file, join(["cell"; string.("n", 1:k)], '\t'))
    for j in 1:length(cellnames)
        print(file, cellnames[j])
        for i in 1:k
            print(file, '\t', database.metadata[neighbors.indexes[i,j]])
        end
        println(file)
    end
end
```

## Installation

First of all, you need to install a Julia compiler.  A recommended way to
install Julia is downloading pre-built binaries. The pre-built binaries of Julia
for several major platforms are distributed at
<https://julialang.org/downloads/>.  Currently, CellFishing.jl supports Julia
0.6 and 0.7 (beta). Then, install CellFishing.jl with the following command:

    $ julia -e 'Pkg.clone("https://github.com/bicycle1885/CellFishing.jl")'

To check the installation, you can try `using CellFishing` in your REPL:

    $ julia
                   _
       _       _ _(_)_     |  A fresh approach to technical computing
      (_)     | (_) (_)    |  Documentation: https://docs.julialang.org
       _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
      | | | | | | |/ _` |  |
      | | |_| | | | (_| |  |  Version 0.7.0-beta.0 (2018-06-24 01:32 UTC)
     _/ |\__'_|_|_|\__'_|  |  Official http://julialang.org/ release
    |__/                   |  x86_64-apple-darwin14.5.0

    julia> using CellFishing  # load the package

    julia>


No error messages mean you have successfully installed CellFishing.jl.


## Command-line interface (WIP)

The bin/cellfishing script is a command-line interface to CellFishing.jl.

    $ ./bin/cellfishing build Plass2018.dge.loom
    Build a search database from Plass2018.dge.loom.
      Loading data ―――――――――――― 13 seconds, 173 milliseconds
      Selecting features ―――――― 1 second, 376 milliseconds
      Creating a database ――――― 16 seconds, 418 milliseconds
      Writing the database ―――― 659 milliseconds
    The serialized database is in Plass2018.dge.loom.cf.
    $ ./bin/cellfishing search Plass2018.dge.loom.cf Plass2018.dge.loom >neighbors.tsv
    Search Plass2018.dge.loom.cf for 10 neighbors.
      Loading the database ―――― 512 milliseconds
      Loading query data ―――――― 12 seconds, 960 milliseconds
      Searching the database ―― 31 seconds, 821 milliseconds
      Writing neighbors ――――――― 64 milliseconds
    $ head -5 neighbors.tsv | cut -f1-3
    plan1_GACTTTCTCTTC      plan1_GACTTTCTCTTC      h2b_TTTTGCTACGGG
    plan1_GTAAGGCGACAN      plan1_GTAAGGCGACAN      gfp_ATTCCTAGCGAT
    plan1_TGGCCCAGCTGC      plan1_TGGCCCAGCTGC      plan1_GACTTTCTCTTC
    plan1_CTCCTGTAATTT      plan1_CTCCTGTAATTT      plan1_ATCCTCCATTAA
    plan1_ATGACGCATAAT      plan1_ATGACGCATAAT      plan1_TACTTGACGGTA

FROM jupyter/datascience-notebook

RUN julia -e 'Pkg.clone("git://github.com/bicycle1885/CellFishing.jl")'

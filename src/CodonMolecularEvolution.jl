module CodonMolecularEvolution

using FASTX, MolecularEvolution, Measures, Compose, PyPlot, StatsBase, Distributions, DataFrames, CSV, NLopt, ParameterHandling, LinearAlgebra

include("shared/shared.jl")
include("difFUBAR/difFUBAR.jl")
include("difFUBAR/grids.jl")
include("../test/benchmark_grid.jl")

# Write your package code here.

end

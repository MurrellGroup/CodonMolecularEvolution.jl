module CodonMolecularEvolution

using FASTX, MolecularEvolution, Measures, Plots, StatsBase, Distributions, DataFrames, CSV, NLopt, ParameterHandling, LinearAlgebra, Phylo

abstract type difFUBARGrid end

include("shared/shared.jl")
include("difFUBAR/difFUBAR.jl")
include("difFUBAR/grids.jl")
include("../test/benchmark_difFUBAR.jl")

export 
    difFUBARBaseline,
    difFUBARParallel,
    difFUBARTreesurgery,
    difFUBARTreesurgeryAndParallel

end

module CodonMolecularEvolution

using FASTX, MolecularEvolution, Measures, Compose, Plots, StatsBase, Distributions, DataFrames, CSV, NLopt, ParameterHandling, LinearAlgebra
using NNlib, Distributions, Zygote, AdvancedHMC, LogDensityProblems, SimpleUnPack, AbstractMCMC, LogDensityProblemsAD

abstract type difFUBARGrid end

include("shared/shared.jl")
include("difFUBAR/difFUBAR.jl")
include("difFUBAR/grids.jl")
include("../test/benchmark_difFUBAR.jl")

include("FUBAR/FUBAR.jl")
include("smoothFUBAR/smoothFUBAR.jl")
include("smoothFUBAR/restrictedFUBAR.jl")
export 
    difFUBARBaseline,
    difFUBARParallel,
    difFUBARTreesurgery,
    difFUBARTreesurgeryAndParallel,
    FUBAR,
    smoothFUBAR,
    restricted_smoothFUBAR

end

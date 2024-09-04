module CodonMolecularEvolution

using FASTX, MolecularEvolution, Measures, Compose, Plots, StatsBase, Distributions, DataFrames, CSV, NLopt, ParameterHandling, LinearAlgebra
using NNlib, Distributions, Zygote, AdvancedHMC, LogDensityProblems, SimpleUnPack, AbstractMCMC, LogDensityProblemsAD, Interpolations, MCMCChains

abstract type difFUBARGrid end
abstract type BAMEgrid end #Bayesian Approaches to Mixed Effects

include("shared/shared.jl")
include("difFUBAR/difFUBAR.jl")
include("difFUBAR/grids.jl")
include("../test/benchmark_difFUBAR.jl")

include("FUBAR/FUBAR.jl")
include("smoothFUBAR/smoothFUBAR.jl")
include("FLAVOR/FLAVOR.jl")

include("MEME/MEME.jl")

include("FAME/FAME.jl")

include("simulations/simulations.jl")

export 
    difFUBARBaseline,
    difFUBARParallel,
    difFUBARTreesurgery,
    difFUBARTreesurgeryAndParallel,
    FUBAR,
    smoothFUBAR

end

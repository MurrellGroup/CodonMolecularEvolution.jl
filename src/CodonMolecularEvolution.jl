module CodonMolecularEvolution

using FASTX, MolecularEvolution, Measures, Compose, Plots, StatsBase, Distributions, DataFrames, CSV, NLopt, ParameterHandling, LinearAlgebra
using NNlib, Distributions, Zygote, AdvancedHMC, LogDensityProblems, SimpleUnPack, AbstractMCMC, LogDensityProblemsAD, DelimitedFiles

abstract type difFUBARGrid end

include("shared/shared.jl")
include("difFUBAR/difFUBAR.jl")
include("difFUBAR/grids.jl")
include("../test/benchmark_difFUBAR.jl")

include("FUBAR/FUBAR.jl")
include("smoothFUBAR/smoothFUBAR.jl")
include("simulations/alphabeta/alphabeta.jl")
export 
    difFUBARBaseline,
    difFUBARParallel,
    difFUBARTreesurgery,
    difFUBARTreesurgeryAndParallel,
    FUBAR,
    smoothFUBAR,
    standard_tree_sim,
    sim_alphabeta_seqs

end

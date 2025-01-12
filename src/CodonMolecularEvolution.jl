module CodonMolecularEvolution

using FASTX, MolecularEvolution, Measures, Compose, Plots, StatsBase, Distributions, DataFrames, CSV, NLopt, ParameterHandling, LinearAlgebra, LaTeXStrings, Random
using NNlib, Distributions, Zygote, AdvancedHMC, LogDensityProblems, SimpleUnPack, AbstractMCMC, LogDensityProblemsAD, Interpolations, MCMCChains

abstract type difFUBARGrid end

include("shared/shared.jl")
include("difFUBAR/difFUBAR.jl")
include("difFUBAR/grids.jl")
include("../test/benchmark_difFUBAR.jl")

include("FUBAR/FUBAR.jl")
include("smoothFUBAR/smoothFUBAR.jl")

include("simulations/alphabeta/alphabeta.jl")
include("simulations/ou_hb.jl")

export 
    difFUBARBaseline,
    difFUBARParallel,
    difFUBARTreesurgery,
    difFUBARTreesurgeryAndParallel,
    FUBAR,
    smoothFUBAR,
    dNdS,
    HBdNdS,
    approx_std2maxdNdS,
    approx_maxdNdS2std,
    HB98AA_matrix,
    ShiftingHBSimModel,
    ShiftingHBSimPartition,
    PiecewiseOUModel,
    shiftingHBviz,
    HBviz
end

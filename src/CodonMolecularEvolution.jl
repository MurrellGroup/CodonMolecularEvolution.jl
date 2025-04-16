module CodonMolecularEvolution

using FASTX, MolecularEvolution, Measures, StatsBase, Distributions, DataFrames, CSV, NLopt, ParameterHandling, LinearAlgebra, Phylo, LaTeXStrings, Random
using NNlib, Distributions, Zygote, AdvancedHMC, LogDensityProblems, SimpleUnPack, AbstractMCMC, LogDensityProblemsAD, Interpolations, MCMCChains
using PDMats, BenchmarkTools, ForwardDiff, Mooncake
using EllipticalSliceSampling
abstract type difFUBARGrid end

include("shared/shared.jl")
include("difFUBAR/difFUBAR.jl")
include("difFUBAR/grids.jl")
include("../test/benchmark_difFUBAR.jl")

include("FUBAR/FUBAR.jl")
include("simulations/alphabeta/alphabeta.jl")
include("simulations/ou_hb.jl")
include("FUBAR/gaussianFUBAR.jl")
include("FUBAR/krylov.jl")
include("FUBAR/grid_utilities.jl")
export 
    difFUBARBaseline,
    difFUBARParallel,
    difFUBARTreesurgery,
    difFUBARTreesurgeryAndParallel,
    FUBAR,
    # smoothFUBAR,
    dNdS,
    HBdNdS,
    std2maxdNdS,
    maxdNdS2std,
    HB98AA_matrix,
    ShiftingHBSimModel,
    ShiftingHBSimPartition,
    PiecewiseOUModel,
    shiftingHBviz,
    HBviz,
    # generate_RJGP_model,
    # reversible_slice_sampling,
    # plot_logposteriors_with_transitions,
    # gridplot,
    # define_gaussian_model,
    # sample_gaussian_model,
    # gaussian_sample_postprocessing,
    FUBAR_analysis,
    SKBDIFUBAR,
    alphabetagrid,
    DirichletFUBAR,
    FIFEFUBAR,
    FUBARgrid
end

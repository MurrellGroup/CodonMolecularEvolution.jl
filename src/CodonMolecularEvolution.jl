module CodonMolecularEvolution

using FASTX, MolecularEvolution, Measures, Plots, StatsBase, Distributions, DataFrames, CSV, NLopt, ParameterHandling, LinearAlgebra, Phylo, LaTeXStrings, Random
using NNlib, Distributions, Zygote, AdvancedHMC, LogDensityProblems, SimpleUnPack, AbstractMCMC, LogDensityProblemsAD, Interpolations, MCMCChains
using PDMats, BenchmarkTools, ForwardDiff, Mooncake
abstract type difFUBARGrid end

include("shared/shared.jl")
include("difFUBAR/difFUBAR.jl")
include("difFUBAR/grids.jl")
include("../test/benchmark_difFUBAR.jl")

include("FUBAR/FUBAR.jl")
include("smoothFUBAR/smoothFUBAR.jl")
include("smoothFUBAR/reversible_slice_sampler.jl")
include("smoothFUBAR/gpFUBAR.jl")
include("smoothFUBAR/matrix_interpolating_gpfubar.jl")
include("smoothFUBAR/wierd_benjamin_trick.jl")
include("simulations/alphabeta/alphabeta.jl")
include("simulations/ou_hb.jl")
include("gaussianFUBAR/gaussianFUBAR.jl")
export 
    difFUBARBaseline,
    difFUBARParallel,
    difFUBARTreesurgery,
    difFUBARTreesurgeryAndParallel,
    FUBAR,
    smoothFUBAR,
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
    generate_RJGP_model,
    reversible_slice_sampling,
    plot_logposteriors_with_transitions,
    gridplot,
    gpFUBAR,
    compute_rjess_to_fubar_permutation,
    non_rj_gpFUBAR,
    kernel_sampling_non_rj_gpFUBAR,
    matrix_interpolating_gp_fubar_HMC_sample,
    ess_benjamin_trick,
    define_gaussian_model,
    sample_gaussian_model
end

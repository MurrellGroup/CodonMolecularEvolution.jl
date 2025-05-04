# This file should include all the stuff that is common across FUBAR methods.
# A FUBAR method at its most basic level should be able to take in a tagged tree
# and output an empirical prior distribution corresponding to a discretisation of the
# bivariate distribution over synonymous/non-synonymous substitution
struct FUBARGrid{T}
    grid_values::Vector{T}
    alpha_vec::Vector{T}
    beta_vec::Vector{T}
    alpha_ind_vec::Vector{Int64}
    beta_ind_vec::Vector{Int64}
    cond_lik_matrix::Matrix{T}
    LL_offset::T
    sites::Int64
    grid_function::Function
    LL_matrix::Matrix{T}
end

# All common functions go here and types go here

abstract type FUBARResults end
abstract type FUBARMethod end
abstract type BayesianFUBARMethod <: FUBARMethod end
struct DefaultBayesianFUBARMethod <: BayesianFUBARMethod end

FUBAR_analysis(method::FUBARMethod; analysis_name="fubar_analysis", write=false) = nothing

function FUBAR_tabulate_results(method::FUBARMethod,results::FUBARResults; analysis_name = "", write = false) end

# Specific functions that are shared

function FUBAR_init(treestring; verbosity=1, exports=true, disable_binarize=false, ladderize_tree=false)

    tree = gettreefromnewick(treestring, FelNode, disable_binarize=disable_binarize)
    if ladderize_tree
        MolecularEvolution.ladderize!(tree)
    end
    verbosity > 0 && println("Step 1: Initialization.")
    return tree
end
sites(p::LazyPartition{CodonPartition}) = p.memoryblocks[1].sites
sites(p::CodonPartition) = p.sites

function FUBAR_grid(tree, GTRmat, F3x4_freqs, code; verbosity=1, grid_function=x -> 10^(x / 6.578947368421053 - 1.502) - 0.0423174293933042, num_grid_points=20)
    verbosity > 0 && println("Step 3: Calculating conditional likelihoods.")

    grid_values = grid_function.(1:num_grid_points)

    LL_matrix = zeros(length(grid_values)^2, sites(tree.message[1]))
    alpha_vec = zeros(length(grid_values)^2)
    alpha_ind_vec = zeros(Int64, length(grid_values)^2)
    beta_vec = zeros(length(grid_values)^2)
    beta_ind_vec = zeros(Int64, length(grid_values)^2)

    i = 1
    for (a, alpha) in enumerate(grid_values)
        for (b, beta) in enumerate(grid_values)
            alpha_vec[i], beta_vec[i] = alpha, beta
            alpha_ind_vec[i], beta_ind_vec[i] = a, b
            m = DiagonalizedCTMC(MolecularEvolution.MG94_F3x4(alpha, beta, GTRmat, F3x4_freqs))
            felsenstein!(tree, m)
            combine!(tree.message[1], tree.parent_message[1])
            LL_matrix[i, :] .= MolecularEvolution.site_LLs(tree.message[1])
            i += 1
        end
    end
    maxi_shift = maximum(LL_matrix, dims=1)
    prob_matrix = exp.(LL_matrix .- maxi_shift)
    sum_shift = sum(prob_matrix, dims=1)
    prob_matrix ./= sum_shift
    LLO = sum(maxi_shift .+ log.(sum_shift))
    return FUBARGrid(grid_values, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec, prob_matrix, LLO, size(prob_matrix, 2), grid_function, LL_matrix)
end

#Packaging "everything before the conditional likelihoods"
function alphabetagrid(seqnames::Vector{String}, seqs, treestring::String;
    verbosity=1, code=MolecularEvolution.universal_code, optimize_branch_lengths=false)
    tree = FUBAR_init(treestring, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit_2steps(seqnames, seqs, tree, x -> x, code, verbosity=verbosity, optimize_branch_lengths=optimize_branch_lengths)
    return FUBAR_grid(tree, GTRmat, F3x4_freqs, code, verbosity=verbosity)
end
function init_path(analysis_name)
    splt = splitpath(analysis_name)[1:end-1]
    if length(splt) > 0
        mkpath(joinpath(splt))
    end
end

struct BayesianFUBARResults{T} <: FUBARResults 
    positive_posteriors::Vector{T}
    purifying_posteriors::Vector{T}
    beta_posterior_mean::Vector{T}
    alpha_posterior_mean::Vector{T}
    posterior_alpha::Matrix{T}
    posterior_beta::Matrix{T}
    posterior_mean::Vector{T}
    theta_chain::Union{Nothing, Vector{Vector{T}}}
    kernel_parameters::Union{Nothing, Vector{Vector{T}}}
end

# For some Bayesian methods, we use EM instead of MCMC and in that case do not get a chain. 

function FUBAR_bayesian_postprocessing(θs::Vector{Vector{T}}, grid::FUBARGrid{T}, kernel_parameters::Vector{Vector{T}}) where {T}
    θ = mean(θs)
    results = FUBAR_bayesian_postprocessing(θ, grid)
    return BayesianFUBARResults(results.positive_posteriors, results.purifying_posteriors, results.beta_posterior_mean, results.alpha_posterior_mean, results.posterior_alpha, results.posterior_beta, θ, θs, kernel_parameters)
end

function FUBAR_bayesian_postprocessing(θ::Vector{T}, grid::FUBARGrid{T}) where {T}
    positive_filter = grid.beta_ind_vec .> grid.alpha_ind_vec
    purifying_filter = grid.beta_ind_vec .< grid.alpha_ind_vec
    weighted_matrix = grid.cond_lik_matrix .* θ
    # Renormalisation
    weighted_matrix = weighted_matrix ./ sum(weighted_matrix, dims=1)
    positive_posteriors = sum(weighted_matrix[positive_filter, :], dims=1)[:]
    purifying_posteriors = sum(weighted_matrix[purifying_filter, :], dims=1)[:]
    beta_posterior_mean = sum(weighted_matrix .* grid.beta_vec, dims=1)[:]
    alpha_posterior_mean = sum(weighted_matrix .* grid.alpha_vec, dims=1)[:]

    n_grid = Int(sqrt(length(grid.beta_vec)))  # should be 20
    weighted_sites = reshape(weighted_matrix, n_grid, n_grid, :)
    # Keep as matrices (n_grid × num_sites)
    posterior_alpha = reshape(sum(weighted_sites, dims=2), n_grid, :)  # Sum over beta dimension
    posterior_beta = reshape(sum(weighted_sites, dims=1), n_grid, :)   # Sum over alpha dimension
    
    return BayesianFUBARResults(positive_posteriors, purifying_posteriors,
        beta_posterior_mean, alpha_posterior_mean,
        posterior_alpha, posterior_beta, θ, nothing, nothing)
end


function FUBAR_tabulate_results(method::DefaultBayesianFUBARMethod, results::BayesianFUBARResults, grid::FUBARGrid; analysis_name = "bayesian_analysis", write = false)
    
    df_results = DataFrame(site=1:size(grid.cond_lik_matrix, 2),
        positive_posterior=results.positive_posteriors,
        purifying_posterior=results.purifying_posteriors,
        beta_posterior_mean=results.beta_posterior_mean,
        alpha_pos_mean=results.alpha_posterior_mean)

    if write
        init_path(analysis_name)
        CSV.write(analysis_name * "_results.csv", df_results)
    end
    return df_results
end

# SKBDI - Smooth Kernel Based Density Inference

# BEGIN DIRICHLET FUBAR
function FUBAR_fitEM(con_lik_matrix, iters, conc; verbosity=1)
    verbosity > 0 && println("Step 4: Model fitting.")
    L = size(con_lik_matrix, 1)
    LDAθ = weightEM(con_lik_matrix, ones(L) ./ L, conc=conc, iters=iters)
    return LDAθ
end

struct DirichletFUBAR <: BayesianFUBARMethod end

"""
    FUBAR_analysis(method::DirichletFUBAR, grid::FUBARGrid{T};
                  analysis_name = "dirichlet_fubar_analysis",
                  write = true,
                  posterior_threshold = 0.95,
                  verbosity = 1,
                  code = MolecularEvolution.universal_code,
                  optimize_branch_lengths = false,
                  concentration = 0.5,
                  iterations = 2500,
                  volume_scaling = 1.0) where {T}

Perform a Fast Unconstrained Bayesian AppRoximation (FUBAR) analysis using a Dirichlet process.

# Arguments
- `method::DirichletFUBAR`: Empty struct to dispatch the original FUBAR method 
- `grid::FUBARGrid{T}`: the FUBARGrid to perform inference on

# Keywords
- `analysis_name::String="dirichlet_fubar_analysis"`: File names
- `write::Bool=true`: Whether to write results to files. Will plot if MolecularEvolutionViz is present
- `posterior_threshold::Float64=0.95`: Posterior probability classification threshold for 
- `verbosity::Int=1`: Control level of output messages (0=none, higher values=more details)
- `code=MolecularEvolution.universal_code`: Molecular code to use
- `optimize_branch_lengths::Bool=false`: ?
- `concentration::Float64=0.5`: Concentration parameter for the Dirichlet process
- `iterations::Int=2500`: Number of EM algorithm iterations
- `volume_scaling::Float64=1.0`: ?

# Returns
- A tuple containing:
  - `df_results`: DataFrame with FUBAR analysis results
  - `θ`: Parameter estimates from the EM algorithm

# Description
Takes in a FUBARGrid object and outputs results for sites obtained from the FUBAR method 
"""
function FUBAR_analysis(method::DirichletFUBAR, grid::FUBARGrid{T}; 
    analysis_name = "dirichlet_fubar_analysis",
    write = true,
    posterior_threshold = 0.95, 
    verbosity = 1,
    code = MolecularEvolution.universal_code,
    optimize_branch_lengths = false,
    concentration = 0.5,
    iterations = 2500,
    volume_scaling = 1.0) where {T}
    
    if write
        init_path(analysis_name)
    end
    
    θ = FUBAR_fitEM(grid.cond_lik_matrix, iterations, concentration, 
                verbosity = verbosity)
                
    results = FUBAR_bayesian_postprocessing(θ, grid)
    
    df_results = FUBAR_tabulate_results(DefaultBayesianFUBARMethod(),results, grid, analysis_name = analysis_name, write = write)

    FUBAR_plot_results(PlotsExtDummy(), method, results, grid, analysis_name = analysis_name, write = write)

    return df_results, (θ = θ, )



end
# END DIRICHLET FUBAR


## HERE BEGINS FIFE FUBAR 
function interpolating_LRS(grid)
    itp = interpolate(grid, BSpline(Cubic(Line(OnGrid()))))
    
    #Null model:
    ab = brents_method_minimize(x -> -itp(x,x), 1, 20, identity, 1e-20)
    LL_null = itp(ab,ab)

    #Alt model:
    a,b,LL_alt = alternating_maximize(itp, (1.0,20.0), (1.0,20.0))
    LRS = 2(LL_alt-LL_null)
    p = 1-cdf(Chisq(1), LRS)
    return p, (p_value = p, LRS = LRS, LL_alt = LL_alt, α_alt = a, β_alt = b, LL_null = LL_null, αβ_null = ab, sel = ifelse(p<0.05,ifelse(b>a, "Positive", "Purifying"), ""))
end
function alternating_maximize(f, a_bounds, b_bounds; final_tol = 1e-20, max_iters = 50)
    a = sum(a_bounds)/2
    b = sum(b_bounds)/2
    m = -f(a,b)
    for i in 1:3
        a = brents_method_minimize(x -> -f(x,b), a_bounds[1], a_bounds[2], identity, 1e-5)
        b = brents_method_minimize(x -> -f(a,x), b_bounds[1], b_bounds[2], identity, 1e-5)
    end
    m = -f(a,b)
    next_m = -Inf
    for i in 1:max_iters
        a = brents_method_minimize(x -> -f(x,b), a_bounds[1], a_bounds[2], identity, final_tol)
        b = brents_method_minimize(x -> -f(a,x), b_bounds[1], b_bounds[2], identity, final_tol)
        next_m = -f(a,b)
        if abs(next_m - m) < final_tol
            break;
        end
        m = next_m
    end
    return a,b,-m
end
struct FrequentistFUBARResults{T} <: FUBARResults
    site_p_value::Vector{T} # p value for α ≠ β at site
    site_LRS::Vector{T} # The raw likelihood ratio
    ha_loglikelihood::Vector{T} # The maximised ll of the alt hyp
    fitted_alpha_ha::Vector{T} # ML estimate of alpha under alt hyp
    fitted_beta_ha::Vector{T} # ML estimate of beta under alt hyp
    hzero_loglikelihood::Vector{T}
    fitted_rate_hzero::Vector{T} # ML estimate of transition rate under null hyp
end
struct FIFEFUBAR <: FUBARMethod end

"""
    FUBAR_analysis(method::FIFEFUBAR, grid::FUBARGrid{T};
                analysis_name = "fife_analysis",
                verbosity = 1,
                write = true,
                positive_tail_only = false) where {T}

Perform a FUBAR type analysis using the FIFE (Frequentist Inference For Evolution) approach.

# Arguments
- `method::FIFEFUBAR`: Empty struct to dispatch on
- `grid::FUBARGrid{T}`: Grid containing data to perform inference on

# Keywords
- `analysis_name::String="fife_analysis"`: Name for the analysis output files and directory
- `verbosity::Int=1`: Control level of output messages (0=none, higher values=more details)
- `write::Bool=true`: Whether to write results to files
- `positive_tail_only::Bool=false`: If true, uses a one-tailed test for positive selection only

# Returns
- `results::FrequentistFUBARResults{T}`: A struct containing the frequentist analysis results

# Description
Frequentist method that gives p-values for site-wise alpha/beta tests

This function performs likelihood ratio tests at each site using interpolated likelihood surfaces.
When `positive_tail_only=true`, the p-values are adjusted to reflect a one-tailed test that only 
considers positive selection (β > α) by using a dirac delta/Chi-square mixture
"""
function FUBAR_analysis(method::FIFEFUBAR, grid::FUBARGrid{T}; 
    analysis_name = "fife_analysis", 
    verbosity = 1, 
    write = true,
    positive_tail_only = false) where {T}
    
    LLmatrix = reshape(grid.LL_matrix, length(grid.grid_values), 
                    length(grid.grid_values), :)
                    
    # Note: dim1 = beta, dim2 = alpha, so we transpose going in:
    stats = [interpolating_LRS(LLmatrix[:,:,i]') for i in 1:size(LLmatrix, 3)]

    # If using one-tailed test for positive selection only
    if positive_tail_only
        # For positive selection:
        # - If α > β (not positive selection): p-value = 0.5
        # - If β > α (potential positive selection): p-value = original p-value / 2
        # This comes from the correct null distribution being a 50:50 mixture of a point mass at 0 and a Chi2(1) distribution
        site_p_value = [s[2].β_alt > s[2].α_alt ? s[2].p_value / 2 : 0.5 for s in stats]
    else
        site_p_value = [s[2].p_value for s in stats]
    end
    
    site_LRS = [s[2].LRS for s in stats]
    ha_loglikelihood = [s[2].LL_alt for s in stats]
    
    fitted_alpha_ha = grid.grid_function.([s[2].α_alt for s in stats])
    fitted_beta_ha = grid.grid_function.([s[2].β_alt for s in stats])
    hzero_loglikelihood = [s[2].LL_null for s in stats]
    fitted_rate_hzero = grid.grid_function.([s[2].αβ_null for s in stats])
    
    results = FrequentistFUBARResults{T}(
        site_p_value,
        site_LRS,
        ha_loglikelihood,
        fitted_alpha_ha,
        fitted_beta_ha,
        hzero_loglikelihood,
        fitted_rate_hzero
    )
    if write
        # save_fubar_results(results, analysis_name = analysis_name)
    end
    
    return results
end



function FUBAR_tabulate_results(method::FIFEFUBAR,results::FrequentistFUBARResults; analysis_name = "fife_analysis", write = false)
    n_sites = length(results.site_p_value)
    
    df_results = DataFrame(
        site = 1:n_sites,
        p_value = results.site_p_value,
        LRS = results.site_LRS,
        LL_alt = results.ha_loglikelihood,
        α_alt = results.fitted_alpha_ha,
        β_alt = results.fitted_beta_ha,
        LL_null = results.hzero_loglikelihood,
        αβ_null = results.fitted_rate_hzero,
        sel = [p < 0.05 ? (b > a ? "Positive" : "Purifying") : "" 
            for (p, a, b) in zip(results.site_p_value, results.fitted_alpha_ha, results.fitted_beta_ha)]
    )

    if write
        init_path(analysis_name)
        CSV.write(analysis_name * "_results.csv", df_results)
    end

    return df_results
end


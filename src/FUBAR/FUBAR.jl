# This file should include all the stuff that is common across FUBAR methods.
# A FUBAR method at its most basic level should be able to take in a tagged tree
# and output an empirical prior distribution corresponding to a discretisation of the
# bivariate distribution over synonymous/non-synonymous substitution
struct FUBARgrid{T}
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

function gridplot(grid::FUBARgrid, θ; title="")
    return nothing
end

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
    return FUBARgrid(grid_values, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec, prob_matrix, LLO, size(prob_matrix, 2), grid_function, LL_matrix)
end

struct FUBARPosterior{T}
    positive_posteriors::Vector{T}
    purifying_posteriors::Vector{T}
    beta_posterior_mean::Vector{T}
    alpha_posterior_mean::Vector{T}
    posterior_alpha::Matrix{T}
    posterior_beta::Matrix{T}
    posterior_mean::Vector{T}
end

function FUBAR_shared_postprocessing(θ, grid::FUBARgrid)
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
    
    return FUBARPosterior(positive_posteriors, purifying_posteriors,
        beta_posterior_mean, alpha_posterior_mean,
        posterior_alpha, posterior_beta, θ)
end
function violin_plot_sites(grid::FUBARgrid, posteriors, posterior_alpha::Matrix,
    posterior_beta::Matrix;
    volume_scaling=1.0,
    posterior_threshold=0.95)
    return nothing
end

function plot_FUBAR_posterior(grid::FUBARgrid, posterior::FUBARPosterior;
    posterior_threshold=0.95,
    volume_scaling=1.0)
    return nothing, nothing, nothing
end

function FUBAR_posterior_to_df(grid::FUBARgrid, posterior::FUBARPosterior;
    analysis_name="fubar_analysis",
    write=false)
    df_results = DataFrame(site=1:size(grid.cond_lik_matrix, 2),
        positive_posterior=posterior.positive_posteriors,
        purifying_posterior=posterior.purifying_posteriors,
        beta_posterior_mean=posterior.beta_posterior_mean,
        alpha_pos_mean=posterior.alpha_posterior_mean)

    if write
        CSV.write(analysis_name * "_results.csv", df_results)
    end
    return df_results
end

function FUBAR_analysis(grid::FUBARgrid, θ;
    analysis_name="fubar_analysis",
    save=false,
    posterior_threshold=0.95,
    volume_scaling=1.0,
    verbosity = 1)
    
    posterior = FUBAR_shared_postprocessing(θ, grid)
    positive_violin_plots, purifying_violin_plots,
    posterior_mean_plot = plot_FUBAR_posterior(grid,
        posterior,
        posterior_threshold=
        posterior_threshold,
        volume_scaling=
        volume_scaling)

    df_results = FUBAR_posterior_to_df(grid, posterior, analysis_name=analysis_name, write=save)

    # Create plots_to_save tuple
    plots_to_save = (
        positive_violin = (positive_violin_plots, analysis_name * "_violin_positive.pdf"),
        purifying_violin = (purifying_violin_plots, analysis_name * "_violin_purifying.pdf"),
        posterior_mean = (posterior_mean_plot, analysis_name * "_posterior_mean.pdf")
    )

    # This will be a no-op if Plots.jl is not loaded
    maybe_save_plots(plots_to_save, save)

    return positive_violin_plots, purifying_violin_plots, posterior_mean_plot, df_results
end

# Default implementation that does nothing
function maybe_save_plots(plots_to_save::NamedTuple, save::Bool)
    return nothing
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

abstract type FUBARMethod end

FUBAR_analysis(method::FUBARMethod; analysis_name="fubar_analysis", save=false) = nothing
# SKBDI - Smooth Kernel Based Density Inference

## HERE BEGINS OLD FUBAR
struct DirichletFUBAR{T} <: FUBARMethod
    grid::FUBARgrid{T}
end

function FUBAR_fitEM(con_lik_matrix, iters, conc; verbosity=1)
    verbosity > 0 && println("Step 4: Model fitting.")
    L = size(con_lik_matrix, 1)
    LDAθ = weightEM(con_lik_matrix, ones(L) ./ L, conc=conc, iters=iters)
    return LDAθ
end
function FUBAR_analysis(method::DirichletFUBAR; analysis_name=
                                "dirichlet_fubar_analysis",
                                save=true,
                                posterior_threshold=0.95, verbosity=1,
                                code=MolecularEvolution.universal_code,
                                optimize_branch_lengths=false,
                                em_parameters=(concentration=0.5, 
                                iterations=2500),
                                volume_scaling=1.0)
    if save
        init_path(analysis_name)
    end
    θ = FUBAR_fitEM(method.grid.cond_lik_matrix, em_parameters.iterations, 
                        em_parameters.concentration, verbosity=verbosity)
    analysis = FUBAR_analysis(method.grid, θ; posterior_threshold =
                                                        posterior_threshold, 
                                                        volume_scaling = 
                                                        volume_scaling,
                                                        save = save,
                                                        verbosity = verbosity)
    return analysis, (θ = θ, )
end
## HERE ENDS OLD FUBAR

## HERE BEGINS FIFE FUBAR 

struct FIFEFUBAR{T} <: FUBARMethod
    grid::FUBARgrid{T}
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

#Generalize this to work with any grid dimensions!
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

#Frequentist interpolating fixed-effects FUBAR
#This needs to have exports like regular FUBAR - plots etc
#We could plot LL surfaces for each site under selection, like this:
#=
    itp = interpolate((LLmatrix[:,:,190]'), BSpline(Cubic(Line(OnGrid()))));
    x = range(1, 20, length=200)
    y = range(1, 20, length=200)
    z = @. itp(x', y)
    contour(x, y, z, levels=30, color=:turbo, fill=true, linewidth = 0, colorbar = false, size = (400,400))
=#
function FUBAR_analysis(method::FIFEFUBAR; analysis_name = "fife_analysis", verbosity=1, save=true)
    save && init_path(analysis_name)
    f = method.grid
    LLmatrix = reshape(f.LL_matrix, length(f.grid_values),length(f.grid_values),:) 
    #Note: dim1 = beta, dim2 = alpha, so we transpose going in:
    stats = [interpolating_LRS(LLmatrix[:,:,i]') for i in 1:size(LLmatrix, 3)]
    df_results = DataFrame([s[2] for s in stats])
    df_results.site = 1:size(LLmatrix, 3)
    df_results.α_alt .= f.grid_function.(df_results.α_alt)
    df_results.β_alt .= f.grid_function.(df_results.β_alt)
    df_results.αβ_null .= f.grid_function.(df_results.αβ_null)
    save && CSV.write(analysis_name * "_results.csv", df_results)
    return df_results
end

## HERE ENDS FIFE FUBAR

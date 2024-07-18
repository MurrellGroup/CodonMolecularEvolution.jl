struct LogPosteriorRestrictedRFF{T} # Quite sure this is not needed since we already have it in smoothFUBAR.jl
    FF::Array{T,2}
    condlik::Array{T,2}
    dim::Int
    unflatten::Function
    posmask::Vector{T}

    # parameters for normal prior on ψ
    μψ::Real
    σψ::Real

    function LogPosteriorRestrictedRFF(gridcoords::Array{T,2}, condlik::Array{T,2}, K::Int, σ::T, unflatten, posmask, μψ::Real, σψ::Real) where T
        W = randn(T, 2, K ÷ 2) * σ * T(2π)
        WtX = W' * gridcoords
        FF = [cos.(WtX); sin.(WtX)]
        new{T}(FF, condlik, K, unflatten, posmask, μψ, σψ)
    end
 end


function restricted_thetas(ω,ψ, LP::LogPosteriorRestrictedRFF)
    

    unrestricted_thetas = softmax((ω' * LP.FF)[:])
    
    interpolated_ψ = softmask2(ψ)

    masked_thetas = (LP.posmask .* interpolated_ψ .* unrestricted_thetas) .+ ((1 .- LP.posmask) .* unrestricted_thetas)



    return masked_thetas ./ sum(masked_thetas)

end

function restricted_LL(ω,ψ, LP::LogPosteriorRestrictedRFF)
    theta_vec = restricted_thetas(ω,ψ, LP)
    return sum(log.(theta_vec' * LP.condlik))
end



function (problem::LogPosteriorRestrictedRFF)(θ)
    ω, ψ = problem.unflatten(θ)
    loglikelihood = restricted_LL(ω,ψ, problem)
    #loglikelihood = 1.0 #if you just want to sample from the prior
    logprior = sum(logpdf.(Normal(0.0,1),ω)) + logpdf(Normal(problem.μψ,problem.σψ),ψ)
    loglikelihood + logprior
end

LogDensityProblems.logdensity(p::LogPosteriorRestrictedRFF, θ) = p(θ)
LogDensityProblems.dimension(p::LogPosteriorRestrictedRFF) = p.dim
LogDensityProblems.capabilities(::Type{LogPosteriorRestrictedRFF}) = LogDensityProblems.LogDensityOrder{0}()

#Helper function to viz theta surface
function thetas(ℓπ::LogPosteriorRestrictedRFF, flatθ)
    ω, ψ = ℓπ.unflatten(flatθ)
    return restricted_thetas(ω,ψ, ℓπ)
end


function softmask2(x::T) where T
    if x < 0
        return T(0)
    else
        return T(log(1+x^3))
    end
end

function restricted_FUBAR_HMCfitRFF(con_lik_matrix, alpha_ind_vec, beta_ind_vec, beta_vec, LL_offset, analysis_name; HMC_samples = 500, K = 50, sigma = 0.03,μψ = 0.5, σψ = 0.5, verbosity=1)
    verbosity > 0 && println("Step 4: Estimating posterior surface by HMC.")

    gridcoords = Matrix{Float64}(hcat(alpha_ind_vec, beta_ind_vec)')
    grid_values = beta_vec[1:Int(sqrt(length(beta_vec)))]
    
    ωeights = randn(K)
    initial_ψ = randn()
    initial_θ = (ω = ωeights, ψ = initial_ψ)
    flat_initial_θ, unflatten = value_flatten(initial_θ)

    ℓπ = LogPosteriorRestrictedRFF(gridcoords, con_lik_matrix, K + 1, sigma, unflatten, beta_ind_vec .> alpha_ind_vec, μψ, σψ)
    num_params = length(flat_initial_θ)
    model = AdvancedHMC.LogDensityModel(LogDensityProblemsAD.ADgradient(Val(:Zygote), ℓπ))
    n_samples, n_adapts = HMC_samples, 200
    δ = 0.8
    sampler = NUTS(δ)
    AdvancedHMC.ProgressMeter.ijulia_behavior(:clear)
    samples = AbstractMCMC.sample(
        model,
        sampler,
        n_adapts + n_samples;
        nadapts = n_adapts,
        initial_params = flat_initial_θ,
    )
     

    ψ_samples = [samples[i].z.θ[K + 1] for i in 1:length(samples)]

    ψ_positive_given_data = sum(ψ_samples .> 0) / length(ψ_samples)
    ψ_negative_given_data = 1 - ψ_positive_given_data
    ψ_positive_prior = 1 - cdf(Normal(μψ, σψ),0)
    ψ_negative_prior = 1 - ψ_positive_prior

    bayes_factor = (ψ_positive_given_data / ψ_negative_given_data ) * (ψ_negative_prior / ψ_positive_prior)

    to_write = "BF="*string(bayes_factor)*"\nP(ψ > 0 | data) = "*string(ψ_positive_given_data)

    write(analysis_name*"_psi_chain_samples.txt",to_write)

    

    posterior_mean_θ = mean([thetas(ℓπ, samples[i].z.θ) for i in 101:length(samples)]);

    # Plotting the MCMC trace of θ values (induced by the mixture of gaussians) for 20 grid points (from the largest to the smallest)
    

    return posterior_mean_θ, samples, ℓπ
end



function generate_grid_stuff(posterior_mean_θ, samples, seqnames, seqs, treestring,ℓπ, analysis_name, outpath;
    n_adapts = 200,pos_thresh=0.95, verbosity=1, exports=true,
     code=MolecularEvolution.universal_code, optimize_branch_lengths=false, K = 50, sigma = 0.03,μψ = 1/2, σψ = 1/2)


    grid_stuff = FUBAR_init2grid(seqnames, seqs, treestring, outpath,
        pos_thresh=pos_thresh, verbosity=verbosity, exports=exports, code=code, optimize_branch_lengths=optimize_branch_lengths)

        
   

      

 

end

function generate_save_plots(posterior_mean_θ, samples, seqnames, seqs, treestring,ℓπ,grid_stuff, analysis_name, outpath;
    n_adapts = 200,pos_thresh=0.95, verbosity=1, exports=true,
     code=MolecularEvolution.universal_code, optimize_branch_lengths=false, K = 50, sigma = 0.03,μψ = 1/2, σψ = 1/2)


    con_lik_matrix, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec, LL_offset = grid_stuff

        gridcoords = Matrix{Float64}(hcat(alpha_ind_vec, beta_ind_vec)')
        grid_values = beta_vec[1:Int(sqrt(length(beta_vec)))]

   

      

    anim = @animate for i ∈ (n_adapts+1):1:length(samples)
        gridplot(alpha_ind_vec,beta_ind_vec,grid_values,thetas(ℓπ, samples[i].z.θ))
    end
    gif(anim, analysis_name * "_posterior_θ_samples.gif", fps = 15)

    lpt_plot = plot(LL_offset .+ [samples[i].stat.log_density for i in (n_adapts+1):length(samples)], label = "Log Posterior", size = (700,350), xlabel = "Iterations")
    savefig(lpt_plot, analysis_name * "_log_posterior_trace.pdf")

    lpt_plot_burnin = plot(LL_offset .+ [samples[i].stat.log_density for i in 1:length(samples)], label = "Log Posterior", size = (700,350), xlabel = "Iterations")
    savefig(lpt_plot_burnin,analysis_name * "_log_posterior_trace_with_burnin.pdf")

    inds_to_plot = sortperm(posterior_mean_θ, rev = true)[(1:(Int(sqrt(length(alpha_ind_vec))))).^2]

    

     posterior_θ_plot = plot(stack([thetas(ℓπ, samples[s].z.θ) for s in (n_adapts+1):length(samples)])[inds_to_plot,:]',
        legend = :none, alpha = 0.5, size = (700,350), xlabel = "Iterations")
     savefig(posterior_θ_plot,analysis_name * "_posterior_θ_trace.pdf")

     ψ_samples = [samples[i].z.θ[K + 1] for i in 1:length(samples)]


     ψ_plot = plot(softmask2.(ψ_samples), label = "Trace of mixing parameter ψ", size = (700,350), xlabel = "Iterations")
     savefig(ψ_plot, analysis_name * "_mixing_parameter_trace.pdf")
 

end

function test_positive_selection_smoothFUBAR(seqnames, seqs, treestring, outpath;
    pos_thresh=0.95, verbosity=1, exports=true, code=MolecularEvolution.universal_code, optimize_branch_lengths=false,
    K = 50, sigma = 0.03, HMC_samples = 500)
   
    con_lik_matrix, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec, LL_offset = FUBAR_init2grid(seqnames, seqs, treestring, outpath,
        pos_thresh=pos_thresh, verbosity=verbosity, exports=exports, code=code, optimize_branch_lengths=optimize_branch_lengths)

    θ, samples, ℓπ = restricted_FUBAR_HMCfitRFF(con_lik_matrix, alpha_ind_vec, beta_ind_vec, beta_vec, LL_offset, outpath, verbosity=verbosity, HMC_samples = HMC_samples, K = K, sigma = sigma)

    verbosity > 0 && println("Step 5: Tabulating results and saving plots.")
    df_results = FUBAR_tabulate_from_θ(con_lik_matrix, θ, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec, outpath, posterior_threshold = pos_thresh)
    #Return df, (tuple of partial calculations needed to re-run tablulate)
    return df_results, (con_lik_matrix, θ, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec), samples, ℓπ
end
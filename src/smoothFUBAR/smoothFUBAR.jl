struct LogPosteriorRFF{T} 
    gridcoords::Array{T,2}
    W::Array{T,2}
    condlik::Array{T,2}
    dim::Int
    unflatten::Function
    posmask::Vector{T}
    function LogPosteriorRFF(gridcoords::Array{T,2}, condlik::Array{T,2}, K::Int, σ::T, unflatten, posmask) where T
        W = randn(T, 2, K ÷ 2) * σ * T(2π)
        new{T}(gridcoords, W, condlik, K, unflatten, posmask)
    end
end

function thetas(ω, LP::LogPosteriorRFF)
    WtX = LP.W'LP.gridcoords
    FF = [cos.(WtX); sin.(WtX)]
    return softmax((ω' * FF)[:])
end

function LL(ω, LP::LogPosteriorRFF)
    theta_vec = thetas(ω, LP)
    return sum(log.(theta_vec' * LP.condlik))
end

#=
#Mixture LL with and without positive selection. MCMC chain doesn't mix nicely.
function LL(ω, LP::LogPosteriorRFF)
    theta_vec = thetas(ω, LP)
    nopos_theta_vec = theta_vec .* nonposmask
    nopos_theta_vec = nopos_theta_vec ./ sum(nopos_theta_vec)
    return logsumexp(sum(log.(theta_vec' * LP.condlik))-0.6931471805599453,sum(log.(nopos_theta_vec' * LP.condlik))-0.6931471805599453)
end
=#

function (problem::LogPosteriorRFF)(θ)
    ω = problem.unflatten(θ)
    loglikelihood = LL(ω, problem)
    #loglikelihood = 1.0 #if you just want to sample from the prior
    logprior = sum(logpdf.(Normal(0.0,1.0),ω))
    loglikelihood + logprior
end

LogDensityProblems.logdensity(p::LogPosteriorRFF, θ) = p(θ)
LogDensityProblems.dimension(p::LogPosteriorRFF) = p.dim
LogDensityProblems.capabilities(::Type{LogPosteriorRFF}) = LogDensityProblems.LogDensityOrder{0}()

#Helper function to viz theta surface
function thetas(ℓπ::LogPosteriorRFF, flatθ)
    ω = ℓπ.unflatten(flatθ)
    return thetas(ω, ℓπ)
end

#Activation function for continuous supression of positive selection
function softmask(x::T) where T
    if x < 0
        return T(0)
    elseif x > 1
        return T(1)
    else
        return T(1/2 + sin(x*pi - pi/2)/2)
    end
end

function FUBAR_HMCfitRFF(con_lik_matrix, alpha_ind_vec, beta_ind_vec, beta_vec, LL_offset, analysis_name; HMC_samples = 500, K = 50, sigma = 0.03, verbosity=1)
    verbosity > 0 && println("Step 4: Estimating posterior surface by HMC.")

    gridcoords = Matrix{Float64}(hcat(alpha_ind_vec, beta_ind_vec)')
    grid_values = beta_vec[1:Int(sqrt(length(beta_vec)))]
    
    ωeights = randn(K)
    initial_θ = (ω = ωeights)
    flat_initial_θ, unflatten = value_flatten(initial_θ)
    ℓπ = LogPosteriorRFF(gridcoords, con_lik_matrix, K, sigma, unflatten, beta_ind_vec .> alpha_ind_vec)
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
    # anim = @animate for i ∈ (n_adapts+1):1:length(samples)
    #     gridplot(alpha_ind_vec,beta_ind_vec,grid_values,thetas(ℓπ, samples[i].z.θ))
    # end
    # gif(anim, analysis_name * "_posterior_θ_samples.mp4", fps = 15)

    # plot(LL_offset .+ [samples[i].stat.log_density for i in (n_adapts+1):length(samples)], label = "Log Posterior", size = (700,350), xlabel = "Iterations")
    # savefig(analysis_name * "_log_posterior_trace.pdf")

    # plot(LL_offset .+ [samples[i].stat.log_density for i in 1:length(samples)], label = "Log Posterior", size = (700,350), xlabel = "Iterations")
    # savefig(analysis_name * "_log_posterior_trace_with_burnin.pdf")

    posterior_mean_θ = mean([thetas(ℓπ, samples[i].z.θ) for i in 101:length(samples)]);

    #Plotting the MCMC trace of θ values (induced by the mixture of gaussians) for 20 grid points (from the largest to the smallest)
    # inds_to_plot = sortperm(posterior_mean_θ, rev = true)[(1:(Int(sqrt(length(alpha_ind_vec))))).^2]
    # plot(stack([thetas(ℓπ, samples[s].z.θ) for s in (n_adapts+1):length(samples)])[inds_to_plot,:]',
    #    legend = :none, alpha = 0.5, size = (700,350), xlabel = "Iterations")
    # savefig(analysis_name * "_posterior_θ_trace.pdf")

    return posterior_mean_θ
end

function smoothFUBAR(seqnames, seqs, treestring, outpath;
    pos_thresh=0.95, verbosity=0, exports=true, code=MolecularEvolution.universal_code, optimize_branch_lengths=false,
    K = 50, sigma = 0.03, HMC_samples = 500)
   
    con_lik_matrix, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec, LL_offset = FUBAR_init2grid(seqnames, seqs, treestring, outpath,
        pos_thresh=pos_thresh, verbosity=verbosity, exports=exports, code=code, optimize_branch_lengths=optimize_branch_lengths)

    θ = FUBAR_HMCfitRFF(con_lik_matrix, alpha_ind_vec, beta_ind_vec, beta_vec, LL_offset, outpath, verbosity=verbosity, HMC_samples = HMC_samples, K = K, sigma = sigma)
    
    verbosity > 0 && println("Step 5: Tabulating results and saving plots.")
    # df_results = FUBAR_tabulate_from_θ(con_lik_matrix, θ, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec, outpath, posterior_threshold = pos_thresh)
    #Return df, (tuple of partial calculations needed to re-run tablulate)
    return Nothing, (con_lik_matrix, θ, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec)
end


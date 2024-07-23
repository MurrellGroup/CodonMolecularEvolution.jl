function plot_trace(samples, param, LP, outpath; transform = x->x, title = "")
    plot([transform(LP.unflatten(s.z.θ)[param]) for s in samples], label = "$(param)", size = (700,350), xlabel = "Iterations", title = title)
    savefig(outpath)
end

#Activation function for continuous supression of positive selection
function sin_softmask(x::T) where T
    if x < 0
        return T(0)
    elseif x > 1
        return T(1)
    else
        return T(1/2 + sin(x*pi - pi/2)/2)
    end
end

function cubic_softmask(x::T) where T
    if x < 0
        return T(0)
    else
        return T(log(1+x^3))
    end
end

abstract type HMC_FUBAR end

struct LogPosteriorRFF{T, H}
    gridcoords::Array{T,2}
    FF::Array{T,2}
    condlik::Array{T,2}
    dim::Int
    unflatten::Function
    flat_initial_θ::Vector{T}
    posmask::Vector{T}
    posneginfmask::Vector{T} #0 where beta>alpha, otherwise -Inf
    nonposneginfmask::Vector{T} #0 where beta<=alpha, otherwise -Inf
    priors::NamedTuple
    function LogPosteriorRFF(method::H, f::FUBARgrid{T}, K::Int, σ::T, initial_θ, priors) where T where H
        @assert Set(keys(priors)) == Set(keys(initial_θ))
        flat_initial_θ, unflatten = value_flatten(initial_θ)
        gridcoords = Matrix{T}(hcat(f.alpha_ind_vec, f.beta_ind_vec)')
        posmask = f.beta_ind_vec .> f.alpha_ind_vec
        W = randn(T, 2, K ÷ 2) * σ * T(2π)
        WtX = W'gridcoords
        FF = [cos.(WtX); sin.(WtX)]
        new{T, H}(gridcoords, FF, f.cond_lik_matrix, length(flat_initial_θ), unflatten, flat_initial_θ, posmask, log.(posmask), log.(1 .- posmask), priors)
    end
end

#Functions shared for all methods:
#Log likelihood
function loglik(θ, LP::LogPosteriorRFF)
    theta_vec = thetas(θ, LP)
    return sum(log.(theta_vec' * LP.condlik))
end

logprior(θ, LP::LogPosteriorRFF) = sum([sum(logpdf.((LP.priors[s],),θ[s])) for s in keys(θ)])

#Helper function to viz theta surface
thetas(ℓπ::LogPosteriorRFF, flatθ) = thetas(ℓπ.unflatten(flatθ), ℓπ)

#Problem setup for AdvancedHMC
function (problem::LogPosteriorRFF)(flatθ)
    θ = problem.unflatten(flatθ)
    return loglik(θ, problem) + logprior(θ, problem)
    #return = logprior(θ, problem) #if you just want to sample from the prior
end
LogDensityProblems.logdensity(p::LogPosteriorRFF, θ) = p(θ)
LogDensityProblems.dimension(p::LogPosteriorRFF) = p.dim
LogDensityProblems.capabilities(::Type{LogPosteriorRFF}) = LogDensityProblems.LogDensityOrder{0}()

#Sampler
function HMCsample(ℓπ, n_samples; n_adapts = 200)
    model = AdvancedHMC.LogDensityModel(LogDensityProblemsAD.ADgradient(Val(:Zygote), ℓπ))
    δ = 0.8
    sampler = NUTS(δ)
    AdvancedHMC.ProgressMeter.ijulia_behavior(:clear)
    samples = AbstractMCMC.sample(
        model,
        sampler,
        n_adapts + n_samples;
        nadapts = n_adapts,
        initial_params = ℓπ.flat_initial_θ,
    )
    return samples
end

#Plots that are commong to all RFF methods
function core_plots(ℓπ, samples, posterior_mean_θ, f, analysis_name, n_adapts)
    anim = @animate for i ∈ (n_adapts+1):1:length(samples)
        gridplot(f.alpha_ind_vec,f.beta_ind_vec,f.grid_values,thetas(ℓπ, samples[i].z.θ))
    end
    gif(anim, analysis_name * "_posterior_θ_samples.gif", fps = 15)

    plot(f.LL_offset .+ [samples[i].stat.log_density for i in (n_adapts+1):length(samples)], label = "Log Posterior", size = (700,350), xlabel = "Iterations")
    savefig(analysis_name * "_log_posterior_trace.pdf")

    plot(f.LL_offset .+ [samples[i].stat.log_density for i in 1:length(samples)], label = "Log Posterior", size = (700,350), xlabel = "Iterations")
    savefig(analysis_name * "_log_posterior_trace_with_burnin.pdf")

    inds_to_plot = sortperm(posterior_mean_θ, rev = true)[(1:(Int(sqrt(length(f.alpha_ind_vec))))).^2]
    plot(stack([thetas(ℓπ, samples[s].z.θ) for s in (n_adapts+1):length(samples)])[inds_to_plot,:]',
        legend = :none, alpha = 0.5, size = (700,350), xlabel = "Iterations")
    savefig(analysis_name * "_posterior_θ_trace.pdf")
end

#Generic other_plots that does nothing - will be overridden when needed
other_plots(LP::LogPosteriorRFF, samples, f, analysis_name, n_adapts) = nothing

#HMC sampling:
function FUBAR_HMCfitRFF(method::HMC_FUBAR, f::FUBARgrid, analysis_name; HMC_samples = 500, n_adapts = 200, K = 50, sigma = 0.03, verbosity=1)
    verbosity > 0 && println("Step 4: Estimating posterior surface by HMC.")
    ℓπ = model_init(method, f, K, sigma)
    samples = HMCsample(ℓπ, HMC_samples, n_adapts = n_adapts)
    posterior_mean_θ = mean([thetas(ℓπ, samples[i].z.θ) for i in n_adapts+1:length(samples)])
    # core_plots(ℓπ, samples, posterior_mean_θ, f, analysis_name, n_adapts)
    # other_plots(ℓπ, samples, f, analysis_name, n_adapts)
    return posterior_mean_θ, samples, ℓπ
end

function pre_init2grid(seqnames, seqs, treestring, outpath;
    pos_thresh=0.95, verbosity=1, exports=true, code=MolecularEvolution.universal_code, optimize_branch_lengths=false)
    f = FUBAR_init2grid(seqnames, seqs, treestring, outpath,
        pos_thresh=pos_thresh, verbosity=verbosity, exports=exports, code=code, optimize_branch_lengths=optimize_branch_lengths)
    
    return f

end

#Main FUBAR call:
function smoothFUBAR(method::HMC_FUBAR, seqnames, seqs, treestring, outpath;
    pos_thresh=0.95, verbosity=1, exports=true, code=MolecularEvolution.universal_code, optimize_branch_lengths=false, K = 150, sigma = 0.03, HMC_samples = 2500)
    f = FUBAR_init2grid(seqnames, seqs, treestring, outpath,
        pos_thresh=pos_thresh, verbosity=verbosity, exports=exports, code=code, optimize_branch_lengths=optimize_branch_lengths)
    θ = FUBAR_HMCfitRFF(method, f, outpath, HMC_samples = HMC_samples, K = K, sigma = sigma, verbosity = verbosity)
    df_results = FUBAR_tabulate_from_θ(θ, f, outpath, posterior_threshold = pos_thresh, verbosity = verbosity)
    return df_results, (θ, f) #(tuple of partial calculations needed to re-run tablulate)
end

function smoothFUBAR_precomputed_f(method::HMC_FUBAR, f, outpath;
    pos_thresh=0.95, verbosity=1, K = 150, sigma = 0.03, HMC_samples = 2500)
    θ, samples, ℓπ = FUBAR_HMCfitRFF(method, f, outpath, HMC_samples = HMC_samples, K = K, sigma = sigma, verbosity = verbosity)
    df_results = FUBAR_tabulate_from_θ(θ, f, outpath, posterior_threshold = pos_thresh, verbosity = verbosity)
    return df_results, (θ, f), samples, ℓπ #(tuple of partial calculations needed to re-run tablulate)
end

#Functions specialized for each method:
#RFF FUBAR - no pos selection control
struct FUBARsmooth <: HMC_FUBAR end
thetas(θ, LP::LogPosteriorRFF{Float64, FUBARsmooth}) = softmax(((θ.ω)' * LP.FF)[:])
model_init(m::FUBARsmooth, f::FUBARgrid{T}, K::Int, σ::T) where T = LogPosteriorRFF(m, f, K, σ, (ω = randn(K),), (ω = Normal(0.0,1.0),))
export FUBARsmooth

#Generics for methods that use a smooth interpolation T(φpre) to control positive selection, where φpre<=0 means no positive selection
abstract type WeightedHMC_FUBAR <: HMC_FUBAR end

function other_plots(LP::LogPosteriorRFF{Float64, M}, samples, f, analysis_name, n_adapts) where M<:WeightedHMC_FUBAR
    prior_possel = 1 - cdf(LP.priors.φpre,0)
    poschain = [LP.unflatten(s.z.θ).φpre > 0 for s in samples]
    pospos = sum(poschain/length(poschain))
    bayesfactor = pospos/(1-pospos) * (1-prior_possel)/prior_possel
    ppstr = "P(φ>0|data)=$(round(pospos,sigdigits=4)). BF = $(round(bayesfactor,sigdigits=4))"
    plot_trace(samples[n_adapts+1:end], :φpre, LP, analysis_name*"_φpre_trace.pdf"; transform = x->x, title = ppstr)
    plot_trace(samples[n_adapts+1:end], :φpre, LP, analysis_name*"_φ_trace.pdf"; transform = transf(M()), title = ppstr)
    println(ppstr)
end
model_init(m::WeightedHMC_FUBAR, f::FUBARgrid{T}, K::Int, σ::T) where T = LogPosteriorRFF(m, f, K, σ, (ω = randn(K), φpre = 0.5),(ω = Normal(0.0,1.0), φpre = Normal(0.5,0.5)))

#Where φ is the proportion of sites where beta>alpha
struct FUBARweightedpos <: WeightedHMC_FUBAR end
transf(M::FUBARweightedpos) = sin_softmask
function thetas(θ, LP::LogPosteriorRFF{Float64, FUBARweightedpos})
    logits = ((θ.ω)' * LP.FF)[:]
    φ = sin_softmask(θ.φpre)
    return φ .* softmax(logits .+ LP.posneginfmask) .+ (1-φ) .* softmax(logits .+ LP.nonposneginfmask)
end
export FUBARweightedpos

#Where φ is multiplied to beta>alpha grids
abstract type PostWeightedHMC_FUBAR <: WeightedHMC_FUBAR end
function thetas(θ, LP::LogPosteriorRFF{Float64, M}) where M<:PostWeightedHMC_FUBAR
    pretheta = softmax(((θ.ω)' * LP.FF)[:])
    T = transf(M())
    φ = T(θ.φpre)
    theta = φ .* pretheta .* LP.posmask .+ pretheta .* (1 .- LP.posmask)
    return theta ./ sum(theta)
end

#...where φ is between 0 and 1
struct FUBARsinacpos <: PostWeightedHMC_FUBAR end
transf(M::FUBARsinacpos) = sin_softmask
export FUBARsinacpos

#...where φ is from 0 to Inf
struct FUBARlogoneplusxcubedpos <: PostWeightedHMC_FUBAR end
transf(M::FUBARlogoneplusxcubedpos) = cubic_softmask
export FUBARlogoneplusxcubedpos



#=
#Mixture LL with and without positive selection. MCMC chain doesn't mix nicely.
function LL(ω, LP::LogPosteriorRFF)
    theta_vec = thetas(ω, LP)
    nopos_theta_vec = theta_vec .* nonposmask
    nopos_theta_vec = nopos_theta_vec ./ sum(nopos_theta_vec)
    return logsumexp(sum(log.(theta_vec' * LP.condlik))-0.6931471805599453,sum(log.(nopos_theta_vec' * LP.condlik))-0.6931471805599453)
end
=#
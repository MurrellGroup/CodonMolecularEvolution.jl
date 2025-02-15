using EllipticalSliceSampling: ESSModel, ESSState, ESS
struct ReversibleSliceSampler
    full_Σ_cholesky::PDMat # Full kernel matrix, Cholesky factor    
    model_dimensions::Vector{Int64} # Dimensions of the models
    model_priors::Vector{Float64} # Priors of the models
    marginal_Σ_cholesky::Vector{PDMat} # Marginal kernel matrices, Cholesky factors
    conditional_cholesky::Dict{Tuple{Int64,Int64},PDMat} # Conditional kernel matrices, Cholesky factors
    loglikelihood::Function # Log likelihood function
    ESS_models::Vector{ESSModel}
end

function generate_reversible_slice_sampler(Σ::Matrix{Float64}, model_dimensions::Vector{Int64}, model_priors::Vector{Float64}, loglikelihood::Function; prior_only=false)
    full_Σ_cholesky = PDMat(cholesky(Symmetric(Σ)))  # Ensure matrix is symmetric
    marginal_Σ_cholesky = [PDMat(cholesky(Symmetric(Σ[1:i, 1:i]))) for i in model_dimensions]
    conditional_cholesky = Dict{Tuple{Int64,Int64},PDMat}()
    for j in 1:length(model_dimensions)
        for i in 1:j
            smaller_model_dimension = model_dimensions[i]
            larger_model_dimension = model_dimensions[j]

            # Extract the conditional covariance matrix properly, this is very unsatisfying since we have 
            # The cholesky factors, but it works and this code only runs once.
            Σ_11 = Σ[1:smaller_model_dimension, 1:smaller_model_dimension]
            Σ_12 = Σ[1:smaller_model_dimension, (smaller_model_dimension+1):larger_model_dimension]
            Σ_22 = Σ[(smaller_model_dimension+1):larger_model_dimension, (smaller_model_dimension+1):larger_model_dimension]

            # Compute conditional covariance
            conditional_cov = Symmetric(Σ_22 - Σ_12' * (Σ_11 \ Σ_12))

            conditional_cholesky[(i, j)] = PDMat(cholesky(conditional_cov))
        end
    end
    ESS_models = [ESSModel(MvNormal(zeros(model_dimensions[i]), marginal_Σ_cholesky[i]), loglikelihood) for i in 1:length(model_dimensions)]

    return ReversibleSliceSampler(full_Σ_cholesky, model_dimensions, model_priors, marginal_Σ_cholesky, conditional_cholesky, prior_only ? x -> 0.0 : loglikelihood, ESS_models)
end

function generate_jump_proposal(sampler::ReversibleSliceSampler, θ::Vector{Float64}, current_model_index::Int64)
    proposed_model_index = rand(setdiff(1:length(sampler.model_dimensions), current_model_index))
    dimension_difference = sampler.model_dimensions[proposed_model_index] - sampler.model_dimensions[current_model_index]
    if dimension_difference > 0
        new_θ = sampler.conditional_cholesky[(current_model_index, proposed_model_index)] * randn(dimension_difference)
        log_jacobian = logdet(sampler.conditional_cholesky[(current_model_index, proposed_model_index)])
        return [θ; new_θ], proposed_model_index, log_jacobian
    else
        new_θ = θ[1:sampler.model_dimensions[proposed_model_index]] # Deterministic projection
        log_jacobian = -logdet(sampler.conditional_cholesky[(proposed_model_index, current_model_index)])
        return new_θ, proposed_model_index, log_jacobian
    end

end

function sample_within_model(sampler::ReversibleSliceSampler, θ::Vector{Float64}, current_model_index::Int64)
    ESS_model = sampler.ESS_models[current_model_index]
    ESS_state = ESSState(θ, sampler.loglikelihood(θ))
    return AbstractMCMC.step(Random.default_rng(), ESS_model, ESS(), ESS_state)[1]
end

function reversible_jump_ess(sampler::ReversibleSliceSampler; n_samples=1000, model_switching_probability=0.05)
    θ = []
    model_indices = []

    current_model_index = rand(1:length(sampler.model_dimensions))
    current_θ = rand(MvNormal(zeros(sampler.model_dimensions[current_model_index]), sampler.marginal_Σ_cholesky[current_model_index]))
    
    ℓπ = (x, model_index) -> sampler.loglikelihood(x) + logpdf(MvNormal(zeros(sampler.model_dimensions[model_index]), sampler.marginal_Σ_cholesky[model_index]), x)[1] + log(sampler.model_priors[model_index])

    logposteriors = [ℓπ(current_θ, current_model_index)]
    while length(θ) < n_samples
        if rand() < model_switching_probability
            proposed_θ, proposed_model_index, log_jacobian = generate_jump_proposal(sampler, current_θ, current_model_index)
            likelihood_ratio = sampler.loglikelihood(proposed_θ) - sampler.loglikelihood(current_θ)

            model_prior_ratio = log(sampler.model_priors[proposed_model_index]) - log(sampler.model_priors[current_model_index])

            # Total log acceptance ratio
            log_α = likelihood_ratio + model_prior_ratio + log_jacobian
            if log(rand()) < log_α
                current_θ = proposed_θ
                current_model_index = proposed_model_index
                push!(θ, current_θ)
                push!(model_indices, current_model_index)
                push!(logposteriors, ℓπ(current_θ, current_model_index))
            end
        else
            current_θ = sample_within_model(sampler, current_θ, current_model_index)
            push!(θ, current_θ)
            push!(model_indices, current_model_index)
            push!(logposteriors, ℓπ(current_θ, current_model_index))
        end
    end
    return θ, model_indices, logposteriors
end
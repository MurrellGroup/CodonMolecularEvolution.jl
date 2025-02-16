using EllipticalSliceSampling: ESSModel, ESSState, ESS

using PDMats
using LinearAlgebra
using Distributions

struct ReversibleSliceSampler
    full_Σ::PDMat # Full kernel matrix, Cholesky factor    
    model_dimensions::Vector{Int64} # Dimensions of the models
    model_priors::Vector{Float64} # Priors of the models
    marginal_Σ::Vector{PDMat} # Marginal kernel matrices, Cholesky factors
    conditional_Σ::Dict{Tuple{Int64,Int64},PDMat} # Conditional kernel matrices, Cholesky factors
    mean_matrices::Dict{Tuple{Int64,Int64},Matrix{Float64}} # Mean matrices
    loglikelihood::Function # Log likelihood function
    ESS_models::Vector{ESSModel}
end

function generate_reversible_slice_sampler(Σ, model_dimensions::Vector{Int64}, model_priors::Vector{Float64}, loglikelihood::Function)
    model_dimensions = model_dimensions
    model_priors = model_priors

    # Convert marginal matrices to PDMat
    marginal_Σ = [PDMat(Σ[1:i, 1:i]) for i in model_dimensions]
    conditional_Σ = Dict{Tuple{Int64,Int64},PDMat}()
    mean_matrices = Dict{Tuple{Int64,Int64},Matrix{Float64}}()
    
    for j in 1:length(model_dimensions)
        for i in 1:(j-1)
            dim_i = model_dimensions[i]
            dim_j = model_dimensions[j]
        
            Σ_11 = Σ[1:dim_i, 1:dim_i] 
            Σ_12 = Σ[1:dim_i, (dim_i+1):dim_j]  
            Σ_21 = Σ[(dim_i+1):dim_j, 1:dim_i]  
            Σ_22 = Σ[(dim_i+1):dim_j, (dim_i+1):dim_j]  

            conditional_Σ[(i, j)] = PDMat(Σ_22 - Σ_21 * inv(Σ_11) * Σ_12)
            mean_matrices[(i, j)] = Σ_21 * inv(Σ_11)
        end
    end
    ESS_models = [ESSModel(MvNormal(zeros(model_dimensions[i]), marginal_Σ_cholesky[i]), loglikelihood) for i in 1:length(model_dimensions)]

    return ReversibleSliceSampler(PDMat(Σ), model_dimensions, model_priors, marginal_Σ, conditional_Σ, mean_matrices, loglikelihood, ESS_models)
end

function generate_jump_proposal(current_model::Int64,current_θ::Vector{Float64}, sampler::ReversibleSliceSampler)
    proposed_model = rand(setdiff(1:length(sampler.model_dimensions), current_model))
    proposed_dimensions = sampler.model_dimensions[proposed_model]
    if proposed_model < current_model
        proposed_θ = current_θ[1:proposed_dimensions]
        log_jacobian = 0.0
        return proposed_θ, proposed_model
    else
        conditional_μ = sampler.mean_matrices[(current_model, proposed_model)] * current_θ
        conditional_Σ = sampler.conditional_Σ[(current_model, proposed_model)]
        new_parameters = rand(MvNormal(conditional_μ, conditional_Σ))
        log_jacobian = 0
        return [current_θ; new_parameters], proposed_model, log_jacobian
    end
end

function reversible_slice_sampling(model::ReversibleSliceSampler; n_samples=1000, n_burnin=200, prior_only = false, jump_proposal_probability = 0.2)
    current_model_index = rand(1:length(model.model_dimensions))
    current_θ = rand(MvNormal(model.model_priors, model.full_Σ))
    samples = []
    model_indicies = []
    logposterios = []
    ℓπ = x -> model.loglikelihood(x) + logpdf(x,MvNormal(zeros(length(x)), model.marginal_Σ[current_model_index]))

    for i in 1:n_samples
       if rand() < jump_proposal_probability
            proposed_θ, proposed_model, log_jacobian = generate_jump_proposal(current_model, current_θ, model)
            likelihood_ratio = model.loglikelihood(proposed_θ) - model.loglikelihood(current_θ)
            prior_ratio = log(model.model_priors[proposed_model]) - log(model.model_priors[current_model])
            # We need the Metropolis-Hastings-Green ratio, which will according to GPT just be likelihood ratio + model prior ratio
            metropolis_hastings_green_ratio = likelihood_ratio + prior_ratio
            if log(rand()) < metropolis_hastings_green_ratio
                current_θ = proposed_θ
                current_model = proposed_model
            end
            push!(samples, current_θ)
            push!(model_indicies, current_model)
            push!(logposterios, ℓπ(current_θ))
       else 
            current_θ = sample_within_model(sampler, current_θ, current_model)
            push!(samples, current_θ)
            push!(model_indicies, current_model)
            push!(logposterios, ℓπ(current_θ))
        end
    end
end

function sample_within_model(sampler::ReversibleSliceSampler, θ::Vector{Float64}, current_model_index::Int64)
    ESS_model = sampler.ESS_models[current_model_index]
    ESS_state = ESSState(θ, sampler.loglikelihood(θ))
    return AbstractMCMC.step(Random.default_rng(), ESS_model, ESS(), ESS_state)[1]
end
using EllipticalSliceSampling: ESSModel, ESSState, ESS
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

function generate_reversible_slice_sampler(Σ, model_dimensions::Vector{Int64}, model_priors::Vector{Float64}, loglikelihood::Function; prior_only = false)
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

            # Add small regularization term to ensure positive definiteness
            conditional_cov = Σ_22 - Σ_21 * inv(Σ_11) * Σ_12
            symmetric_conditional_cov = Symmetric((conditional_cov + conditional_cov')/2)   
            regularized_cov = symmetric_conditional_cov + I * 1e-8  # Add small constant to diagonal
            conditional_Σ[(i, j)] = PDMat(regularized_cov)
            mean_matrices[(i, j)] = Σ_21 * inv(Σ_11)
        end
    end
    ESS_models = [ESSModel(MvNormal(zeros(model_dimensions[i]), marginal_Σ[i]), loglikelihood) for i in 1:length(model_dimensions)]
    actual_loglikelihood = prior_only ? x -> 0 : loglikelihood
    return ReversibleSliceSampler(PDMat(Σ), model_dimensions, model_priors, marginal_Σ, conditional_Σ, mean_matrices, actual_loglikelihood, ESS_models)
end

function generate_jump_proposal(current_model::Int64,current_θ::Vector{Float64}, sampler::ReversibleSliceSampler)
    
    if current_model == length(sampler.model_dimensions)
        proposed_model = current_model - 1
    elseif current_model == 1
        proposed_model = current_model + 1
    else
        proposed_model = current_model + rand([-1, 1])
    end

    proposed_dimensions = sampler.model_dimensions[proposed_model]
    if proposed_model < current_model
        proposed_θ = current_θ[1:proposed_dimensions]
        # This corrects for a non-symmetric proposal distribution
        correction_factor = current_model == length(sampler.model_dimensions) ? log(0.5) : 0.0 
        return proposed_θ, proposed_model, correction_factor
    else
        conditional_μ = sampler.mean_matrices[(current_model, proposed_model)] * current_θ
        conditional_Σ = sampler.conditional_Σ[(current_model, proposed_model)]
        new_parameters = rand(MvNormal(conditional_μ, conditional_Σ))
        # This corrects for a non-symmetric proposal distribution
        correction_factor = current_model == 1 ? log(0.5) : 0.0
        return [current_θ; new_parameters], proposed_model, correction_factor
    end
end

function reversible_slice_sampling(model::ReversibleSliceSampler; n_samples=1000, n_burnin=200, prior_only = false, jump_proposal_probability = 0.2)
    current_model_index = rand(1:length(model.model_dimensions))
    current_θ = rand(MvNormal(zeros(model.model_dimensions[current_model_index]), model.marginal_Σ[current_model_index]))
    samples = []
    model_indicies = []
    logposterios = []
    ℓπ = x -> model.loglikelihood(x) + logpdf(MvNormal(zeros(length(x)), model.marginal_Σ[current_model_index]), x)

    for i in 1:n_samples
       if rand() < jump_proposal_probability
            proposed_θ, proposed_model_index, correction_factor = generate_jump_proposal(current_model_index, current_θ, model)
            likelihood_ratio = model.loglikelihood(proposed_θ) - model.loglikelihood(current_θ)
            prior_ratio = log(model.model_priors[proposed_model_index]) - log(model.model_priors[current_model_index])
            metropolis_hastings_green_ratio = likelihood_ratio + prior_ratio + correction_factor    
            if log(rand()) < metropolis_hastings_green_ratio
                current_θ = proposed_θ
                current_model_index = proposed_model_index
            end
            push!(samples, current_θ)
            push!(model_indicies, current_model_index)
            push!(logposterios, ℓπ(current_θ))
       else 
            current_θ = sample_within_model(model, current_θ, current_model_index)
            push!(samples, current_θ)
            push!(model_indicies, current_model_index)
            push!(logposterios, ℓπ(current_θ))
        end
    end
    return samples, model_indicies, logposterios
end

function sample_within_model(sampler::ReversibleSliceSampler, θ::Vector{Float64}, current_model_index::Int64)
    ESS_model = sampler.ESS_models[current_model_index]
    ESS_state = ESSState(θ, sampler.loglikelihood(θ))
    result = AbstractMCMC.step(Random.default_rng(), ESS_model, ESS(), ESS_state)[1]
    return result
end
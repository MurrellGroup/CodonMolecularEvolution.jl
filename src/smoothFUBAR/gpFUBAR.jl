struct RJGPModel
    grid::FUBARgrid
    Σ::Matrix{Float64} # Full kernel
    dimension::Int64
    purifying_prior::Float64
    invp::Vector{Int64}
end

function create_rjess_permutation(N::Int)
    # First collect indices for lower triangular + diagonal (column by column)
    lower = Int[]
    for j in 1:N
        for i in j:N
            idx = (j-1)*N + i
            push!(lower, idx)
        end
    end
    
    # Then collect indices for upper triangular (column by column)
    upper = Int[]
    for j in 1:N
        for i in 1:(j-1)
            idx = (j-1)*N + i
            push!(upper, idx)
        end
    end
    
    # Combine to get permutation
    p = vcat(lower, upper)
    
    # Compute inverse permutation
    invp = similar(p)
    for i in eachindex(p)
        invp[p[i]] = i
    end
    
    return invp
end

# Scuffed fix but hey this might work?
function loglikelihood(model::RJGPModel, θ)
    full_θ = [softmax(θ); zeros(model.dimension - length(θ))]
    return sum(log.(full_θ[model.invp]'model.grid.cond_lik_matrix)) # Here we hopåefully re-permute to FUBAR format
end

function rearrange_kernel_matrix(Σ)
    N = Int64(sqrt(size(Σ)[1]))
    result = zeros(N*N, N*N)
    
    # Create mapping from FUBAR format to RJESS format
    rjess_to_fubar = Int[]
    
    # First get indices for lower triangular + diagonal elements in FUBAR format
    for j in 1:N
        for i in j:N
            idx = (j-1)*N + i  # FUBAR column-major index
            push!(rjess_to_fubar, idx)
        end
    end
    
    # Then get indices for upper triangular elements in FUBAR format
    for j in 1:N
        for i in 1:(j-1)
            idx = (j-1)*N + i  # FUBAR column-major index
            push!(rjess_to_fubar, idx)
        end
    end
    
    # Rearrange covariance matrix according to this mapping
    for i in 1:length(rjess_to_fubar)
        for j in 1:length(rjess_to_fubar)
            result[i,j] = Σ[rjess_to_fubar[i], rjess_to_fubar[j]]
        end
    end
    
    return result
end
function gaussian_kernel_matrix(grid; kernel_scaling=1.0)
    return kernel_matrix(grid, distance_function=x -> exp(-x / (2 * kernel_scaling^2)))
end

function kernel_matrix(grid; distance_function=x -> exp(-x))
    n_points = length(grid.alpha_ind_vec)
    distances = zeros(n_points, n_points)
    for i in 1:n_points
        for j in 1:n_points
            distance = (grid.alpha_ind_vec[i] - grid.alpha_ind_vec[j])^2 +
                       (grid.beta_ind_vec[i] - grid.beta_ind_vec[j])^2
            distances[i, j] = distance_function(distance)
        end
    end
    return distances
end

function generate_RJGP_model(grid::FUBARgrid; kernel_scaling = 1.0, purifying_prior = 1/2)
    Σ = rearrange_kernel_matrix(gaussian_kernel_matrix(grid, kernel_scaling = kernel_scaling))
    dimension = size(Σ)[1]
    return RJGPModel(grid, Σ, dimension, purifying_prior, create_rjess_permutation(Int64(sqrt(dimension))))
end

function reversible_slice_sampling(model::RJGPModel; ϵ = 0.01, n_samples = 1000, model_switching_probability = 0.3, prior_only = false)
    ll = prior_only ? x -> 0 : x -> loglikelihood(model, x)
    full_covariance = 1/2 * (model.Σ + model.Σ') + (ϵ * I) # Tikhonov + posdef, need big ϵ for high kernel scaling
    N = Int64(sqrt(model.dimension))
    smallest_model_dimension = Int64(N * (N + 1) / 2)
    model_dimensions = collect(smallest_model_dimension:10:(N^2)) # Hard-coded 10 is a bit ugly, we will fix this later
    # Uniform prior for positive selection models, "lump" prior for purifying
    puryfing_model_prior = (1 - model.purifying_prior) / (length(model_dimensions)-1)
    model_priors = [[model.purifying_prior]; puryfing_model_prior .* ones(length(model_dimensions)-1)]
    println("Model dimension: ",model_dimensions[end])
    println("Σ minimum eigenvalue: ",minimum(eigvals(full_covariance)))
    problem = RJESSProblem(ll, full_covariance, model_dimensions, model_priors)

    return rj_ess(problem, n_samples = n_samples, model_switching_probability = model_switching_probability)
end

function plot_logposteriors_with_transitions(model_indices, logposteriors)
    # Initialize empty plot
    p = plot(legend=true)
    
    # Find transition points
    transitions = findall(diff(model_indices) .!= 0)
    
    # Add first segment
    start_idx = 1
    
    for transition in transitions
        # Plot segment up to transition
        plot!(p, start_idx:transition, logposteriors[start_idx:transition], 
              linewidth=2, label=(start_idx==1 ? "Log Posterior" : false))
        
        # Add vertical red line at transition
        vline!([transition], color=:red, 
               linestyle=:solid, label=(start_idx==1 ? "Model Transition" : false))
        
        # Start next segment after transition
        start_idx = transition + 1
    end
    
    # Plot final segment
    if start_idx <= length(logposteriors)
        plot!(p, start_idx:length(logposteriors), 
              logposteriors[start_idx:end], 
              linewidth=2, label=false)
    end
    
    return p
end
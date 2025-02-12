struct RJGPModel
    grid::FUBARgrid
    Σ::Matrix{Float64} # Full kernel
    dimension::Int64
    purifying_prior::Float64
    invp::Vector{Int64}
end

function create_rjess_permutation(N::Int)
    rjess_to_fubar = Int[]

    # Upper triangular in RJESS format first
    for col in 2:N
        for row in 1:(col-1)  # Top down for each column above diagonal
            fubar_idx = (col - 1) * N + row
            push!(rjess_to_fubar, fubar_idx)
        end
    end

    # Lower triangular + diagonal in RJESS format
    for col in 1:N
        for row in N:-1:col  # Bottom to top, including diagonal
            fubar_idx = (col - 1) * N + row
            push!(rjess_to_fubar, fubar_idx)
        end
    end

    # Create inverse permutation
    invp = similar(rjess_to_fubar)
    for i in eachindex(rjess_to_fubar)
        invp[rjess_to_fubar[i]] = i
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
    result = zeros(N * N, N * N)

    # Create mapping from RJESS indices to FUBAR indices
    rjess_to_fubar = Int[]

    # First map upper triangular elements
    for col in 2:N
        for row in 1:(col-1)
            fubar_idx = (col - 1) * N + row
            push!(rjess_to_fubar, fubar_idx)
        end
    end

    # Then map lower triangular + diagonal elements
    # Following exact pattern from RJESS grid
    for col in 1:N
        for row in N:-1:col  # Bottom to top, including diagonal
            fubar_idx = (col - 1) * N + row
            push!(rjess_to_fubar, fubar_idx)
        end
    end

    # Rearrange kernel matrix using this mapping
    for (i, fi) in enumerate(rjess_to_fubar)
        for (j, fj) in enumerate(rjess_to_fubar)
            result[i, j] = Σ[fi, fj]
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

function generate_RJGP_model(grid::FUBARgrid; kernel_scaling=1.0, purifying_prior=1 / 2)
    Σ = rearrange_kernel_matrix(gaussian_kernel_matrix(grid, kernel_scaling=kernel_scaling))
    dimension = size(Σ)[1]
    return RJGPModel(grid, Σ, dimension, purifying_prior, create_rjess_permutation(Int64(sqrt(dimension))))
end

function reversible_slice_sampling(model::RJGPModel; ϵ=0.01, n_samples=1000, model_switching_probability=0.3, prior_only=false)
    ll = prior_only ? x -> 0 : x -> loglikelihood(model, x)
    full_covariance = 1 / 2 * (model.Σ + model.Σ') + (ϵ * I) # Tikhonov + posdef, need big ϵ for high kernel scaling
    N = Int64(sqrt(model.dimension))
    smallest_model_dimension = Int64(N * (N + 1) / 2)
    model_dimensions = collect(smallest_model_dimension:10:(N^2)) # Hard-coded 10 is a bit ugly, we will fix this later
    # Uniform prior for positive selection models, "lump" prior for purifying
    puryfing_model_prior = (1 - model.purifying_prior) / (length(model_dimensions) - 1)
    model_priors = [[model.purifying_prior]; puryfing_model_prior .* ones(length(model_dimensions) - 1)]
    println("Model dimension: ", model_dimensions[end])
    println("Σ minimum eigenvalue: ", minimum(eigvals(full_covariance)))
    problem = RJESSProblem(ll, full_covariance, model_dimensions, model_priors)

    return rj_ess(problem, n_samples=n_samples, model_switching_probability=model_switching_probability)
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
            linewidth=2, label=(start_idx == 1 ? "Log Posterior" : false))

        # Add vertical red line at transition
        vline!([transition], color=:red,
            linestyle=:solid, label=(start_idx == 1 ? "Model Transition" : false))

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

function compute_rjess_to_fubar_permutation(v::AbstractVector{Float64}, N::Int)
    result = zeros(N * N)
    idx = N * (N - 1) ÷ 2 + 1  # Start idx for lower triangle + diagonal

    # First process lower triangle + diagonal going up each column
    for col in 1:N
        for row in N:-1:col  # Bottom to top, including diagonal
            fubar_idx = (col - 1) * N + row
            result[fubar_idx] = v[idx]
            idx += 1
        end
    end

    # Then process ALL upper triangular elements
    idx = 1  # Reset idx to start of vector for upper triangle
    for col in 2:N
        for row in 1:(col-1)  # Top down for each column above diagonal
            fubar_idx = (col - 1) * N + row
            result[fubar_idx] = v[idx]
            idx += 1
        end
    end

    return result
end

function softmax(x)
    exp_x = exp.(x .- maximum(x))
    return exp_x ./ sum(exp_x)
end

function format_sample(sample, max_dimension)
    return [softmax(sample); zeros(max_dimension - length(sample))]
end

function compute_purifying_bayes_factor(model_indices, purifying_prior)
    # Count visits to model 1 (purifying) and models > 1 (non-purifying)
    purifying_visits = count(x -> x == 1, model_indices)
    non_purifying_visits = count(x -> x > 1, model_indices)
    
    # Compute posterior odds
    posterior_odds = purifying_visits / non_purifying_visits
    
    # Compute prior odds
    prior_odds = purifying_prior / (1 - purifying_prior)
    
    # Bayes factor = posterior odds / prior odds
    bayes_factor = posterior_odds / prior_odds
    
    return bayes_factor
end

function gpFUBAR(problem::RJGPModel; ϵ=0.01, n_samples=1000, model_switching_probability=0.01)
    samples, model_indices, logposteriors, jump_history = reversible_slice_sampling(problem, ϵ=ϵ, n_samples=n_samples, model_switching_probability=model_switching_probability)

    bayes_factor = compute_purifying_bayes_factor(model_indices, problem.purifying_prior)
    println("Bayes factor (M1/M>1): ", bayes_factor)

    formatted_samples = [format_sample(sample, problem.dimension) for sample in samples]
    fubar_samples = [compute_rjess_to_fubar_permutation(formatted_sample, Int64(sqrt(problem.dimension))) for formatted_sample in formatted_samples]

    posterior_mean = mean(fubar_samples)
    active_parameters = length.(samples)

    # Create individual plots
    posterior_mean_plot = gridplot(problem.grid.alpha_ind_vec, problem.grid.beta_ind_vec, 
                                 problem.grid.grid_values, posterior_mean, 
                                 title="Posterior Mean")
    active_parameter_trace = plot(active_parameters, 
                                title="Active Parameters",
                                xlabel="Iteration",
                                ylabel="Number of Parameters")
    logposterior_plot = plot_logposteriors_with_transitions(model_indices, logposteriors)
    
    # Create diagnostic plots layout
    diagnostic_plots = plot(active_parameter_trace, logposterior_plot,
                          layout=(2,1),
                          size=(600, 800))
    
    # Create animation
    anim = @animate for i in 1:100:length(formatted_samples)
        gridplot(problem.grid.alpha_ind_vec, problem.grid.beta_ind_vec, problem.grid.grid_values, fubar_samples[i])
    end

    return posterior_mean_plot, diagnostic_plots, anim
end

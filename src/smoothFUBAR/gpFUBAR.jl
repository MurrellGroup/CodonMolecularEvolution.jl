# Fix 1: Modify RJGPModel to store the correct permutation
struct RJGPModel
    grid::FUBARgrid
    Σ::Matrix{Float64} # Full kernel
    dimension::Int64
    purifying_prior::Float64
    ess_to_fubar_perm::Vector{Int64}  # Changed name to be explicit
end

"""
    generate_ess_indices(N::Int)

Generate ESS format indices for an NxN grid. Indices are generated along diagonals
from bottom-right to top-left, with each diagonal numbered from bottom to top.
"""
function generate_ess_indices(N::Int)
    indices = zeros(Int, N, N)
    current_index = 1
    
    # Loop over each diagonal, starting from bottom right
    for diag in (2N-2):-1:0
        # Calculate range for this diagonal
        col_start = max(0, diag - N + 1)
        col_end = min(diag, N-1)
        
        # Store positions for this diagonal (bottom to top)
        positions = [(diag - col, col) 
                    for col in col_end:-1:col_start]
        
        # Fill in indices
        for pos in positions
            indices[pos[1] + 1, pos[2] + 1] = current_index
            current_index += 1
        end
    end
    
    return indices
end

"""
    generate_fubar_indices(N::Int)

Generate FUBAR format indices for an NxN grid. Indices are generated column-wise
from bottom to top, starting from the leftmost column.
"""
function generate_fubar_indices(N::Int)
    indices = zeros(Int, N, N)
    current_index = 1
    
    # Fill column by column, bottom to top
    for col in 1:N
        for row in N:-1:1
            indices[row, col] = current_index
            current_index += 1
        end
    end
    
    return indices
end

function get_ess_to_fubar_permutation(N::Int)
    ess = generate_ess_indices(N)
    fubar = generate_fubar_indices(N)
    
    # Create mapping from ESS indices to FUBAR indices
    n_elements = N * N
    perm = zeros(Int, n_elements)
    
    for i in 1:N
        for j in 1:N
            ess_idx = ess[i,j]
            fubar_idx = fubar[i,j]
            # We want to map FROM ess TO fubar, so:
            perm[fubar_idx] = ess_idx  # This line was wrong before!
        end
    end
    
    return perm
end

function get_fubar_to_ess_permutation(N::Int)
    ess = generate_ess_indices(N)
    fubar = generate_fubar_indices(N)
    
    # Create mapping from FUBAR indices to ESS indices
    n_elements = N * N
    perm = zeros(Int, n_elements)
    
    for i in 1:N
        for j in 1:N
            ess_idx = ess[i,j]
            fubar_idx = fubar[i,j]
            # We want to map FROM fubar TO ess, so:
            perm[ess_idx] = fubar_idx  # This line was wrong before!
        end
    end
    
    return perm
end

# Fix 3: Update the loglikelihood function
function loglikelihood(model::RJGPModel, θ)
    full_θ = [softmax(θ); zeros(model.dimension - length(θ))]
    return sum(log.(full_θ[model.ess_to_fubar_perm]'model.grid.cond_lik_matrix))
end


function rearrange_kernel_matrix(Σ)
    fubar_to_ess = get_fubar_to_ess_permutation(Int64(sqrt(size(Σ)[1])))
    return Σ[fubar_to_ess, fubar_to_ess]
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

# Fix 2: Update the model constructor
function generate_RJGP_model(grid::FUBARgrid; kernel_scaling=1.0, purifying_prior=1/2)
    Σ = rearrange_kernel_matrix(gaussian_kernel_matrix(grid, kernel_scaling=kernel_scaling))
    dimension = size(Σ)[1]
    return RJGPModel(
        grid, 
        Σ, 
        dimension, 
        purifying_prior, 
        get_ess_to_fubar_permutation(Int64(sqrt(dimension)))  # Changed to ESS→FUBAR
    )
end
function rjess(model::RJGPModel; ϵ=0.01, n_samples=1000, model_switching_probability=0.3, prior_only=false, diagnostics = false)
    ll = prior_only ? x -> 0 : x -> loglikelihood(model, x)
    full_covariance = 1 / 2 * (model.Σ + model.Σ') + (ϵ * I) # Tikhonov + posdef, need big ϵ for high kernel scaling
    N = Int64(sqrt(model.dimension))
    smallest_model_dimension = Int64(N * (N - 1) / 2)
    model_dimensions = accumulate((x,i) -> x + (N-i), 0:(N-1); init=smallest_model_dimension)
    # Uniform prior for positive selection models, "lump" prior for purifying
    diversifying_model_prior = (1 - model.purifying_prior) / (length(model_dimensions) - 1)
    model_priors = [[model.purifying_prior]; diversifying_model_prior .* ones(length(model_dimensions) - 1)]
    if diagnostics
        println("Model dimension: ", model_dimensions[end])
        println("Σ minimum eigenvalue: ", minimum(eigvals(full_covariance)))
        println("Model priors: ", model_priors)
        println("Model dimensions: ", model_dimensions)
    end
    problem = generate_reversible_slice_sampler(full_covariance, model_dimensions, model_priors, ll; prior_only=prior_only)

    return reversible_slice_sampling(problem, n_samples=n_samples, jump_proposal_probability=model_switching_probability)
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
    ess_to_fubar = get_ess_to_fubar_permutation(N)  # Get ESS→FUBAR permutation
    return v[ess_to_fubar]  # Apply it to transform ESS→FUBAR
 end
 

function softmax(x)
    exp_x = exp.(x .- maximum(x))
    return exp_x ./ sum(exp_x)
end

function format_sample(sample, max_dimension)
    return [softmax(sample); zeros(max_dimension - length(sample))]
end

function compute_purifying_bayes_factor(model_indices, purifying_prior)
    # Calculate burnin (20% of samples)
    burnin = Int(floor(0.2 * length(model_indices)))
    
    # Use only post-burnin samples
    post_burnin_indices = model_indices[(burnin + 1):end]
    
    # Count visits to model 1 (purifying) and models > 1 (non-purifying)
    purifying_visits = count(x -> x == 1, post_burnin_indices)
    non_purifying_visits = count(x -> x > 1, post_burnin_indices)
    
    # Compute posterior odds
    posterior_odds = purifying_visits / non_purifying_visits
    
    # Compute prior odds
    prior_odds = purifying_prior / (1 - purifying_prior)
    
    # Bayes factor = posterior odds / prior odds
    bayes_factor = posterior_odds / prior_odds
    
    return bayes_factor
end

function gpFUBAR(problem::RJGPModel; ϵ=0.01, n_samples=1000, model_switching_probability=0.01, prior_only=false, diagnostics = false)
    samples, model_indices, logposteriors = rjess(problem, ϵ=ϵ, n_samples=n_samples, model_switching_probability=model_switching_probability, prior_only=prior_only)

    # Use all samples for Bayes factor calculation and model frequencies
    bayes_factor = compute_purifying_bayes_factor(model_indices, problem.purifying_prior)
    
    # Calculate and print model time distributions using all samples
    model_counts = countmap(model_indices)
    total_samples = length(model_indices)
    # println("Jump acceptance rate: ", jump_diagnostics.accepted / jump_diagnostics.proposed)
    if diagnostics 
        println("\nModel time distributions:")
    end
        for (model, count) in sort(collect(model_counts))
        percentage = (count / total_samples) * 100
       if diagnostics
        println("Model $model: $(round(percentage, digits=2))% ($(count) samples)")
       end
    end
    println("\nBayes factor (M1/M>1): ", bayes_factor)

    # formatted_samples = [format_sample(sample, problem.dimension) for sample in samples]
    # fubar_samples = [compute_rjess_to_fubar_permutation(formatted_sample, Int64(sqrt(problem.dimension))) for formatted_sample in formatted_samples]

    # posterior_mean = mean(fubar_samples[1:100:end]) # Only plot every 1000th sample
    # active_parameters = length.(samples[1:100:end]) # Only plot every 1000th sample

    # Create individual plots
   #  posterior_mean_plot = gridplot(problem.grid.alpha_ind_vec, problem.grid.beta_ind_vec, 
  #                                problem.grid.grid_values, posterior_mean, 
  #                                title="Posterior Mean")
 #    active_parameter_trace = plot(active_parameters, 
  #                               title="Active Parameters",
  #                               xlabel="Iteration",
  #                               ylabel="Number of Parameters")
# logposterior_plot = plot_logposteriors_with_transitions(model_indices, logposteriors)
    
    # Create diagnostic plots layout
   # diagnostic_plots = plot(active_parameter_trace, logposterior_plot,
   #                       layout=(2,1),
   #                       size=(600, 800))
    
    # Create animation
    # anim = @animate for i in 1:200:length(formatted_samples)
    #     gridplot(problem.grid.alpha_ind_vec, problem.grid.beta_ind_vec, problem.grid.grid_values, fubar_samples[i])
    # end
end

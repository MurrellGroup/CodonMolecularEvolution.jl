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
    println("Model indices: ", model_indices[1:20])
    println("Logposteriors: ", logposteriors[1:20])
    println("Samples: ", samples[1])
    # Use all samples for Bayes factor calculation and model frequencies
    bayes_factor = 1 / compute_purifying_bayes_factor(model_indices, problem.purifying_prior)
    
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

    formatted_samples = [format_sample(sample, problem.dimension) for sample in samples]
    fubar_samples = [compute_rjess_to_fubar_permutation(formatted_sample, Int64(sqrt(problem.dimension))) for formatted_sample in formatted_samples]

    posterior_mean = mean(fubar_samples[1:100:end]) # Only plot every 1000th sample
    active_parameters = length.(samples[1:100:end]) # Only plot every 1000th sample

    # Create individual plots
    posterior_mean_plot = gridplot(problem.grid.alpha_ind_vec, problem.grid.beta_ind_vec, 
                                 problem.grid.grid_values, posterior_mean, 
                                 title="Posterior Mean")
    active_parameter_trace = plot(active_parameters, 
                                 title="Active Parameters",
                                 xlabel="Iteration",
                                 ylabel="Number of Parameters")
    logposterior_plot = plot_logposteriors_with_transitions(model_indices[1:100:end], logposteriors[1:100:end])
    
    # Create diagnostic plots layout with modified dimensions
    diagnostic_plots = plot(active_parameter_trace, logposterior_plot,
                           layout=(2,1),
                           size=(1000, 800))  # Changed from (600, 1800) to (800, 400)
    
    # Create animation
    anim = @animate for i in 1:1000:length(formatted_samples)
        gridplot(problem.grid.alpha_ind_vec, problem.grid.beta_ind_vec, problem.grid.grid_values, fubar_samples[i])
    end
    return posterior_mean_plot, diagnostic_plots, anim
end


function quintic_smooth_transition(x, alpha, beta)
    if x <= alpha
        return 0.0
    elseif x >= beta
        return 1.0
    else
        # Normalize x to [0,1] range
        t = (x - alpha) / (beta - alpha)
        # Quintic polynomial that has zero derivatives at t=0 and t=1
        # and exactly reaches 0 at t=0 and 1 at t=1
        return t * t * t * (10.0 + t * (-15.0 + 6.0 * t))
    end
end

# Dimensions gives a vector of the models incremental dimensionality
function supress_vector(θ, αs, βs, dimensions, smoothing_function)
    grid_θ = θ[1:(end-1)]
    supression_parameter = θ[end]

    softmax_grid_θ = softmax(grid_θ)


    for i in 2:length(dimensions)
        start_index = dimensions[i-1] + 1
        end_index = dimensions[i]
        softmax_grid_θ[start_index:end_index] = softmax_grid_θ[start_index:end_index] .* smoothing_function(supression_parameter, αs[i], βs[i])
    end
    return softmax_grid_θ ./ sum(softmax_grid_θ)
end




function supression_loglikelihood(model::RJGPModel, θ, αs, βs, dimensions)
    
    softmax_grid_θ = supress_vector(θ, αs, βs, dimensions, quintic_smooth_transition)
    return sum(log.(softmax_grid_θ[model.ess_to_fubar_perm]'model.grid.cond_lik_matrix))
end

function non_rj_gpFUBAR(problem::RJGPModel; suppression_σ = 1, n_samples=1000, prior_only=false, diagnostics = false)
    new_Σ = zeros(problem.dimension + 1, problem.dimension + 1)
    new_Σ[1:problem.dimension, 1:problem.dimension] = problem.Σ
    new_Σ[problem.dimension + 1, problem.dimension + 1] = suppression_σ^2
    dimensions = accumulate((x,i) -> x + (20-i), 0:(20-1); init=190) #Hard coded 20,190 for now
    αs = 0.1 .* [i for i in 0:(length(dimensions)-1)]
    βs = 0.1 .* [i for i in 1:(length(dimensions))]

    actual_loglikelihood = prior_only ? x -> 0 : x -> supression_loglikelihood(problem, x, αs, βs, dimensions)

    ESS_model = ESSModel(MvNormal(zeros(problem.dimension + 1), new_Σ), actual_loglikelihood)
    samples = AbstractMCMC.sample(ESS_model, ESS(),n_samples,progress=true)

    

    grid_samples = [compute_rjess_to_fubar_permutation(supress_vector(samples[i], αs, βs, dimensions, quintic_smooth_transition), Int64(sqrt(problem.dimension))) for i in eachindex(samples)]
    
    supression_parameter_samples = [samples[i][end] for i in eachindex(samples)]
    posterior_mean = mean(grid_samples[1:100:end]) # Only plot every 1000th sample
    
    # Create individual plots
    posterior_mean_plot = gridplot(problem.grid.alpha_ind_vec, problem.grid.beta_ind_vec, 
                                 problem.grid.grid_values, posterior_mean, 
                                 title="Posterior Mean")
    supression_parameter_trace = plot(supression_parameter_samples[1:100:end], 
                                 title="Supression Parameter",
                                 xlabel="Iteration",
                                 ylabel="Supression Parameter")
    println("Bayes factor: ", sum(supression_parameter_samples .> 0.0)/sum(supression_parameter_samples .< 0.0))

    anim = @animate for i in 1:1000

        sample = compute_rjess_to_fubar_permutation(supress_vector(samples[i], αs, βs, dimensions, quintic_smooth_transition), Int64(sqrt(problem.dimension)))

        gridplot(problem.grid.alpha_ind_vec, problem.grid.beta_ind_vec, problem.grid.grid_values, sample, title = "Supression Parameter: $(i/500)")
    end
    
    
    return posterior_mean_plot, supression_parameter_trace, anim
end

function kernel_sampling_non_rj_gpFUBAR(problem::RJGPModel; suppression_σ = 1, n_samples=1000, prior_only=false, diagnostics = false, sample_c = true, c_step_size = 0.01, c_initial = 2.0)
    # Initialize c with the provided value
    c = c_initial
    
    # Create initial kernel matrix
    Σ = zeros(problem.dimension + 1, problem.dimension + 1)
    
    # Define dimensions, αs, and βs (same as in non_rj_gpFUBAR)
    dimensions = accumulate((x,i) -> x + (20-i), 0:(20-1); init=190) #Hard coded 20,190 for now
    αs = 0.1 .* [i for i in 0:(length(dimensions)-1)]
    βs = 0.1 .* [i for i in 1:(length(dimensions))]
    
    # Create the initial kernel matrix
    if sample_c
        # If we're sampling c, create the kernel based on current c value
        raw_kernel = kernel_matrix(problem.grid, distance_function = x -> exp(-x / (2*c^2)))
        Σ[1:problem.dimension, 1:problem.dimension] = rearrange_kernel_matrix(raw_kernel)
    else
        # If not sampling c, use the kernel from the problem (like in non_rj_gpFUBAR)
        Σ[1:problem.dimension, 1:problem.dimension] = problem.Σ
    end
    Σ[problem.dimension + 1, problem.dimension + 1] = suppression_σ^2
    
    # Add small regularization to ensure numerical stability
    Σ[1:problem.dimension, 1:problem.dimension] += 1e-6 * I
    
    # Define likelihood function
    actual_loglikelihood = prior_only ? x -> 0 : x -> supression_loglikelihood(problem, x, αs, βs, dimensions)
    
    # Initialize θ with a draw from the prior
    θ = rand(MvNormal(zeros(problem.dimension + 1), Σ))
    
    # Create chains to store results
    θs = [θ]
    cs = [c]
    
    # Track logposteriors
    θ_logposteriors = Float64[]
    c_logposteriors = Float64[]
    
    # Track Metropolis acceptance for c
    metropolis_proposed = 0
    metropolis_accepted = 0
    
    # Calculate initial logposterior for θ
    θ_logprior = logpdf(MvNormal(zeros(problem.dimension + 1), Σ), θ)
    θ_loglikelihood = actual_loglikelihood(θ)
    push!(θ_logposteriors, θ_logprior + θ_loglikelihood)
    
    # Calculate initial logposterior for c
    c_logprior = logpdf(Normal(2, 0.25), c)
    push!(c_logposteriors, c_logprior)
    
    # Set up initial ESS state
    ESS_model = ESSModel(MvNormal(zeros(problem.dimension + 1), Σ), actual_loglikelihood)
    ESS_state = ESSState(θ, θ_loglikelihood)
    
    # Run MCMC
    for i in 1:n_samples
        if diagnostics && i % 100 == 0
            println("Iteration $i of $n_samples")
        end
        
        # First, sample c if requested
        if sample_c
            metropolis_proposed += 1
            c_proposal = c + c_step_size * randn()
            
            # Create new kernel matrix with proposed c
            Σ_proposal = zeros(problem.dimension + 1, problem.dimension + 1)
            raw_kernel_proposal = kernel_matrix(problem.grid, distance_function = x -> exp(-x / (2*c_proposal^2)))
            Σ_proposal[1:problem.dimension, 1:problem.dimension] = rearrange_kernel_matrix(raw_kernel_proposal)
            Σ_proposal[problem.dimension + 1, problem.dimension + 1] = suppression_σ^2
            Σ_proposal[1:problem.dimension, 1:problem.dimension] += 1e-6 * I  # Add regularization
            
            # This is the transport map from our current Σ to the proposed Σ

            A = sqrt(Σ) * inv(sqrt(sqrt(Σ) * Σ_proposal * sqrt(Σ))) * sqrt(Σ)

            proposed_θ = A * θ

            # Calculate Metropolis ratio
            new_θ_logprior = logpdf(MvNormal(zeros(problem.dimension + 1), Σ_proposal), proposed_θ)
            proposed_c_logprior = logpdf(Normal(2, 0.25), c_proposal)
            current_θ_logprior = logpdf(MvNormal(zeros(problem.dimension + 1), Σ), θ)
            current_c_logprior = logpdf(Normal(2, 0.25), c)

            new_θ_loglikelihood = actual_loglikelihood(proposed_θ)
            current_θ_loglikelihood = actual_loglikelihood(θ)

            likelihood_ratio = new_θ_loglikelihood - current_θ_loglikelihood
            metropolis_ratio = proposed_c_logprior - current_c_logprior + new_θ_logprior - current_θ_logprior + likelihood_ratio
            
            if log(rand()) < metropolis_ratio
                metropolis_accepted += 1
                c = c_proposal
                Σ = Σ_proposal
                θ = proposed_θ
                # Update ESS model with new kernel
                ESS_model = ESSModel(MvNormal(zeros(problem.dimension + 1), Σ), actual_loglikelihood)
                
                # Update c logposterior
                push!(c_logposteriors, proposed_c_logprior)
                push!(θ_logposteriors, new_θ_logprior + new_θ_loglikelihood)
                push!(θs, θ)
            else
                # Keep the same c but still record its logposterior
                push!(c_logposteriors, current_c_logprior)
                push!(θ_logposteriors, current_θ_logprior + current_θ_loglikelihood)
                push!(θs, θ)
            end
        end
        
        push!(cs, c)
        
        # Now perform a single ESS step
        result, new_state = AbstractMCMC.step(Random.default_rng(), ESS_model, ESS(), ESS_state)
        θ = result
        ESS_state = new_state
        
        push!(θs, θ)
        
        # Calculate and store new logposterior for θ
        θ_logprior = logpdf(MvNormal(zeros(problem.dimension + 1), Σ), θ)
        θ_loglikelihood = new_state.loglikelihood
        push!(θ_logposteriors, θ_logprior + θ_loglikelihood)
    end
    
    # Print Metropolis acceptance ratio if c was sampled
    if sample_c
        acceptance_ratio = metropolis_accepted / metropolis_proposed
        println("Metropolis acceptance ratio for c parameter: $(round(acceptance_ratio * 100, digits=2))% ($(metropolis_accepted)/$(metropolis_proposed))")
    end
    
    # Process results for visualization and analysis
    grid_samples = [compute_rjess_to_fubar_permutation(
        supress_vector(θs[i], αs, βs, dimensions, quintic_smooth_transition), 
        Int64(sqrt(problem.dimension))
    ) for i in eachindex(θs)]
    
    # Calculate posterior mean using every 100th sample
    posterior_mean = mean(grid_samples[1:100:end])
    
    # Create plots
    posterior_mean_plot = gridplot(
        problem.grid.alpha_ind_vec, 
        problem.grid.beta_ind_vec, 
        problem.grid.grid_values, 
        posterior_mean, 
        title="Posterior Mean"
    )
    
    c_trace_plot = plot(
        cs[1:100:end], 
        title="Kernel Parameter (c)",
        xlabel="Iteration",
        ylabel="c value"
    )
    
    logposterior_plot = plot(
        θ_logposteriors[1:100:end],
        title="Log Posterior Trace",
        xlabel="Iteration",
        ylabel="Log Posterior"
    )
    
    # Create animation
    anim = @animate for i in 1:100:min(1000, length(θs))
        sample = grid_samples[i]
        
        gridplot(
            problem.grid.alpha_ind_vec, 
            problem.grid.beta_ind_vec, 
            problem.grid.grid_values, 
            sample, 
            title="Sample at iteration $i, c = $(round(cs[i], digits=3))"
        )
    end
    
    # Return results
    return (
        posterior_mean_plot = posterior_mean_plot,
        c_trace_plot = c_trace_plot,
        logposterior_plot = logposterior_plot,
        animation = anim,
        θ_chain = (samples = θs, logposteriors = θ_logposteriors),
        c_chain = (samples = cs, logposteriors = c_logposteriors),
        posterior_mean = posterior_mean,
        metropolis_stats = sample_c ? (proposed = metropolis_proposed, accepted = metropolis_accepted, ratio = metropolis_accepted/metropolis_proposed) : nothing
    )
end
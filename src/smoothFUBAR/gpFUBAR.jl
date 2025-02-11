struct RJGPModel
    grid::FUBARgrid
    Σ::Matrix{Float64} # Full kernel
    dimension::Int64
    purifying_prior::Float64
end


function loglikelihood(model::RJGPModel, θ)
    full_θ = [zeros(model.dimension - length(θ)); θ]
    return sum(log.(softmax(full_θ)'model.grid.cond_lik_matrix))
end

function rearrange_kernel_matrix(Σ)
    N = Int64(sqrt(size(Σ)[1]))
    upper = [(i - 1) * N + j for i in 1:N for j in 1:N if i < j]
    lower = [(i - 1) * N + j for i in 1:N for j in 1:N if i >= j]
    p = vcat(upper, lower)
    return Σ[p, p]
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
    return RJGPModel(grid, Σ, dimension, purifying_prior)
end

function reversible_slice_sampling(model::RJGPModel; ϵ = 0.01, n_samples = 1000, model_switching_probability = 0.3, prior_only = false)
    ll = prior_only ? x -> 0 : x -> loglikelihood(model, x)
    full_covariance = 1/2 * (model.Σ + model.Σ') + (ϵ * I) # Tikhonov + posdef, need big ϵ for high kernel scaling
    N = Int64(sqrt(model.dimension))
    smallest_model_dimension = Int64(N * (N + 1) / 2)
    model_dimensions = collect(smallest_model_dimension:(N^2))
    # Uniform prior for positive selection models, "lump" prior for purifying
    puryfing_model_prior = (1 - model.purifying_prior) / (length(model_dimensions)-1)
    model_priors = [[model.purifying_prior]; puryfing_model_prior .* ones(length(model_dimensions)-1)]
    println("Model dimension: ",model_dimensions[end])
    println("Σ minimum eigenvalue: ",minimum(eigvals(full_covariance)))
    problem = RJESSProblem(ll, full_covariance, model_dimensions, model_priors)

    return rj_ess(problem, n_samples = n_samples, model_switching_probability = model_switching_probability)
end

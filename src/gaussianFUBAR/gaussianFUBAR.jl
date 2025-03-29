include("grid_utilities.jl")
include("krylov.jl")

struct AmbientESSProblem{T}
    loglikelihood::Function
    distance_function::Function
    kernel_function::Function # has signature (d, kernelparams)
    gaussian_dimension::Int64 # Actual problem dimension
    kernel_parameter_dimension::Int64 # How many parameters to sample in the kernel
end

function generate_distance_matrix(dimension::Int64, 
                                    distance_function::Function)
    distance_matrix = zeros(dimension, dimension)
    for i in 1:dimension
        for j in 1:i
            d = distance_function(i,j)
            distance_matrix[i, j] = d
            distance_matrix[j, i] = d
        end
    end
    return distance_matrix
end

function generate_distance_matrix(problem::AmbientESSProblem)
    return generate_distance_matrix(problem.gaussian_dimension, 
                                    problem.distance_function)
end

# Samples in ambient space. m is krylov subspace dimension
function transform_sample(problem::AmbientESSProblem, θ::AbstractVector; 
                                                        m = 10)
    distance_matrix = generate_distance_matrix(problem)
    return transform_sample(problem, θ, distance_matrix, m=m)
end
# Epsilon is Tykhonoff regularisation
function transform_sample(problem::AmbientESSProblem, θ::AbstractVector, 
                                        distance_matrix::AbstractMatrix; m=10, ϵ = 1e-6)
    # All kernel parameters are at the end
    kernel_parameters = θ[(problem.gaussian_dimension+1):end] 
    kernel_function = x -> problem.kernel_function(x, kernel_parameters...)
    kernel_matrix = kernel_function.(distance_matrix) + ϵ * I
    # println("Check: The kernel matrix is symmetric: ",issymmetric(kernel_matrix))
    # println("Check: The kernel matrix is positive definite: ",isposdef(kernel_matrix))
    # println("Spectrum of the kernel matrix: ",eigvals(kernel_matrix))
    ambient_gaussian_parameters = θ[1:problem.gaussian_dimension]
    transformed_θ = krylov_sqrt_times_vector(kernel_matrix, 
                                            ambient_gaussian_parameters, m=m)                                       
    return transformed_θ
end

# This function takes in a kernel function and a likelihood function
# Samples in ambient/uncorrelated space, transforms the sample using the resulting
# kernel matrix and performs Elliptical Slice sampling
# It uses Krylov subspace methods 

function kernel_sampling_ess(problem::AmbientESSProblem; m=10, 
                                                        prior_only = false,
                                                        n_samples = 1000, 
                                                        burnin = 200,
                                                        progress = false)
    distance_matrix = generate_distance_matrix(problem)
    loglikelihood = prior_only ? θ -> 0 : θ -> problem.loglikelihood(
                                                        transform_sample(
                                                        problem, θ, 
                                                        distance_matrix, m=m))
    # We have to have μ0=0 since otherwise we would have to re-define the prior 
    # for every choice of kernel parameter
    total_dimension = problem.gaussian_dimension + 
                                problem.kernel_parameter_dimension
    μ0 = zeros(total_dimension)
    Σ0 = diagm(ones(total_dimension))
    prior = MvNormal(μ0, Σ0)
    ESS_model = ESSModel(prior, loglikelihood)
    
    ambient_samples = AbstractMCMC.sample(ESS_model, ESS(), n_samples, 
                                                    progress = progress)

    sample_transformation_function = θ -> transform_sample(problem, θ, 
                                                            distance_matrix, 
                                                                        m=m)
    post_burnin_samples = ambient_samples[(burnin+1):end]
    transformed_samples = sample_transformation_function.(post_burnin_samples)
    return transformed_samples
end

abstract type AmbientProblemModel end

struct SupressionType{T}
    alphas::AbstractVector{T}
    betas::AbstractVector{T}
    dimensions::AbstractVector{Int64}
    transition_function::Function
end

struct GaussianFUBARModel{T} <: AmbientProblemModel
    grid::FUBARgrid{T}
    distance_function::Function
    kernel_function::Function
    gaussian_dimension::Int64
    kernel_parameter_dimension::Int64
    # Permutes the grid structure in case the LL function requires it
    fubar_to_ambient_permutation_vector::Vector{Int64}
    ambient_to_fubar_permutation_vector::Vector{Int64}
    supression_type::SupressionType{T}
end

function define_ambient_problem(model::AmbientProblemModel) end

function supress_vector(supression_type::SupressionType, θ::Vector{Float64})
    T = eltype(θ)
    # The last index is the supression parameter, always
    grid_θ = @view θ[1:(end-1)] 
    y = θ[end]

    s = softmax(grid_θ)
    M = ones(T, length(s))

    @inbounds for i in 2:length(supression_type.dimensions)
        start_index = supression_type.dimensions[i-1] + 1
        end_index = supression_type.dimensions[i]
        f_i = supression_type.transition_function(y, supression_type.alphas[i],
                                                    supression_type.betas[i])
        M[start_index:end_index] .= f_i
    end

    A = s .* M
    T_sum = sum(A)
    return A ./ T_sum
end

function gaussian_fubar_loglikelihood(model::GaussianFUBARModel, 
                                            θ::Vector{Float64})
    grid_θ = supress_vector(model.supression_type, θ)
    permuted_grid_θ = grid_θ[model.ambient_to_fubar_permutation_vector]
    return sum(log.(permuted_grid_θ'model.grid.cond_lik_matrix))
end

function define_ambient_problem(model::GaussianFUBARModel)
    distance_matrix = generate_distance_matrix(model.gaussian_dimension, 
                                                model.distance_function)

    fubar_to_ambient_permutation = vcat(model.fubar_to_ambient_permutation_vector, 
                                                [model.gaussian_dimension])

    permuted_distance_matrix = distance_matrix[
                                fubar_to_ambient_permutation,
                                fubar_to_ambient_permutation]


    ambient_distance_function = (i,j) -> permuted_distance_matrix[i,j]
    
    loglikelihood = θ -> gaussian_fubar_loglikelihood(model, θ)
    return AmbientESSProblem{Float64}(loglikelihood, ambient_distance_function, 
                                            model.kernel_function,
                                            model.gaussian_dimension,
                                            model.kernel_parameter_dimension)
end

function sample_gaussian_model(model::GaussianFUBARModel; m=10, 
                                                        prior_only = false,
                                                        n_samples = 1000, 
                                                        burnin = 200,
                                                        progress = false)
    ambient_problem = define_ambient_problem(model)
    return kernel_sampling_ess(ambient_problem, m=m, prior_only=prior_only, 
                                n_samples=n_samples, burnin=burnin,
                                progress=progress)
end

function standard_fubar_distance_function(grid::FUBARgrid, i,j)
    # This accounts for the scaling parameter
    if (i > length(grid.alpha_ind_vec) || j > length(grid.alpha_ind_vec))
        return i == j ? 0 : Inf
    end

    return (grid.alpha_ind_vec[i] - grid.alpha_ind_vec[j])^2 + 
            (grid.beta_ind_vec[i] - grid.beta_ind_vec[j])^2
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

function define_gaussian_model(grid::FUBARgrid; 
                        distance_function = standard_fubar_distance_function,
                        kernel_function = (d,c) -> exp(-d/(c^2)),
                        kernel_parameter_dimension = 1,
                        supression_type = nothing)

    passed_distance_function = (i,j) -> distance_function(grid, i,j)
    # The part of the model that does not impact the covariance structure
    # is (all of the gridpoints) + (the supression parameter)
    gaussian_dimension = length(grid.alpha_ind_vec) + 1

    grid_dimension = length(grid.grid_values)

    fubar_to_ambient = get_fubar_to_ambient_permutation(grid_dimension)
    ambient_to_fubar = get_ambient_to_fubar_permutation(grid_dimension)

    if isnothing(supression_type)
        last_lower_triangular_index = grid_dimension * (grid_dimension - 1) / 2
        supression_dimensions = accumulate((x, i) -> x + (grid_dimension - i),
                                                        0:(grid_dimension-1); 
                                            init=last_lower_triangular_index)       
        alphas = 0.1 .* [i for i in 0:(length(supression_dimensions)-1)]
        betas = 0.1 .* [i for i in 1:(length(supression_dimensions))]

        supression_type = SupressionType{Float64}(alphas, betas, 
                                                    supression_dimensions,
                                                    quintic_smooth_transition)
    end

    return GaussianFUBARModel{Float64}(grid, passed_distance_function, 
                                    kernel_function,
                                    gaussian_dimension, 
                                    kernel_parameter_dimension, 
                                    fubar_to_ambient, ambient_to_fubar,
                                    supression_type)
end
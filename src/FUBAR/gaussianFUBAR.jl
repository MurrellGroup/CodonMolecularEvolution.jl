struct AmbientESSProblem
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
            d = distance_function(i, j)
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
    m=10)
    distance_matrix = generate_distance_matrix(problem)
    return transform_sample(problem, θ, distance_matrix, m=m)
end
# Epsilon is Tykhonoff regularisation

function krylov_sqrt_times_vector(A,v; m = 15)
    lanczos_iterator = LanczosIterator(A,v)
    factorization = initialize(lanczos_iterator)
    while length(factorization) < m
        expand!(lanczos_iterator, factorization)
    end
    Qm = basis(factorization) 
    Qm_matrix = hcat([Qm[i] for i in 1:length(Qm)]...)
    Tm = rayleighquotient(factorization)  
    Tm_sqrt_firstcol = sqrt(Tm)[:,1]
    return norm(v) * (Qm_matrix * Tm_sqrt_firstcol)
end


function transform_sample(problem::AmbientESSProblem, θ::AbstractVector,
    distance_matrix::AbstractMatrix; m=15, ϵ=1e-6)
    # All kernel parameters are at the end
    kernel_parameters = θ[(problem.gaussian_dimension+1):end]
    kernel_function = x -> problem.kernel_function(x, kernel_parameters...)
    kernel_matrix = kernel_function.(distance_matrix) + ϵ * I
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
    prior_only=false,
    n_samples=1000,
    burnin=200,
    progress=false)

    distance_matrix = generate_distance_matrix(problem)
    loglikelihood = prior_only ? θ -> 0 : θ -> problem.loglikelihood(
        transform_sample(
            problem, θ,
            distance_matrix, m=m))
    # We have to have μ0=0 since otherwise we would have to re-define the prior 
    # for every choice of kernel parameter, 
    total_dimension = problem.gaussian_dimension +
                    problem.kernel_parameter_dimension
    μ0 = zeros(total_dimension)
    Σ0 = diagm(ones(total_dimension))
    prior = MvNormal(μ0, Σ0)
    ESS_model = ESSModel(prior, loglikelihood)

    ambient_samples = AbstractMCMC.sample(ESS_model, ESS(), n_samples,
        progress=progress)

    sample_transformation_function = θ -> transform_sample(problem, θ,
        distance_matrix,
        m=m)
    post_burnin_samples = ambient_samples[(burnin+1):end]
    transformed_samples = sample_transformation_function.(post_burnin_samples)
    kernel_parameter_samples = [ambient_samples[i][(problem.gaussian_dimension+1):end] for i in eachindex(ambient_samples)]
    return transformed_samples, kernel_parameter_samples
end

abstract type AmbientProblemModel end

struct SupressionType{T}
    alphas::AbstractVector{T}
    betas::AbstractVector{T}
    dimensions::AbstractVector{Int64}
    transition_function::Function
end

struct GaussianFUBARModel{T} <: AmbientProblemModel
    grid::FUBARGrid{T}
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


    ambient_distance_function = (i, j) -> permuted_distance_matrix[i, j]

    loglikelihood = θ -> gaussian_fubar_loglikelihood(model, θ)
    return AmbientESSProblem(loglikelihood, ambient_distance_function,
        model.kernel_function,
        model.gaussian_dimension,
        model.kernel_parameter_dimension)
end



function standard_fubar_distance_function(grid::FUBARGrid, i, j)
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

function define_gaussian_model(grid::FUBARGrid;
    distance_function=standard_fubar_distance_function,
    kernel_function=(d, c) -> exp(-d / (c^2)),
    kernel_parameter_dimension=1,
    supression_type=nothing)

    passed_distance_function = (i, j) -> distance_function(grid, i, j)
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

function sample_gaussian_model(model::GaussianFUBARModel; m=10,
    prior_only=false,
    n_samples=1000,
    burnin=200,
    progress=false)

    ambient_problem = define_ambient_problem(model)
    return kernel_sampling_ess(ambient_problem, m=m, prior_only=prior_only,
        n_samples=n_samples, burnin=burnin,
        progress=progress)
end

function gaussian_sample_postprocessing(model::GaussianFUBARModel, θs; thinning=100, m=10)
    thinned_samples = θs[1:thinning:end]
    grid_samples = [supress_vector(model.supression_type, θ)[model.ambient_to_fubar_permutation_vector] for θ in thinned_samples]
    return grid_samples
end

## HERE BEGINS INTEGRATION WIH THE FUBAR INTERACE
struct SKBDIFUBAR <: BayesianFUBARMethod end

"""
    FUBAR_analysis(method::SKBDIFUBAR, grid::FUBARGrid{T};
                analysis_name = "skbdi_fubar_analysis",
                volume_scaling = 1.0,
                write = true,
                verbosity = 1,
                posterior_threshold = 0.95,
                distance_function = standard_fubar_distance_function,
                kernel_function = (d, c) -> exp(-d / c^2),
                kernel_parameter_dimension = 1,
                supression_type = nothing,
                m = 15,
                ϵ = 1e-6,
                n_samples = 1000,
                burnin = 200,
                thinning = 50) where {T}

Perform a Fast Unconstrained Bayesian AppRoximation (FUBAR) analysis using the SKBDI (Smooth Kernel Bayesian Density Inference) approach.

# Arguments
- `method::SKBDIFUBAR`: Empty struct used for dispatch
- `grid::FUBARGrid{T}`: Grid to perform inference on

# Keywords
- `analysis_name::String="skbdi_fubar_analysis"`: Name for the analysis output files and directory
- `volume_scaling::Float64=1.0`: ?
- `write::Bool=true`: Whether to write results to files
- `verbosity::Int=1`: Control level of output messages (0=none, higher values=more details)
- `posterior_threshold::Float64=0.95`: Posterior probability threshold for classification
- `distance_function=standard_fubar_distance_function`: Function used to calculate distances between grid points
- `kernel_function=(d, c) -> exp(-d / c^2)`: Kernel function used for the covariance matrix. 
- `kernel_parameter_dimension::Int=1`: How many kernel parameters are taken in by the kernel bandwidth function. 
- `supression_type=nothing`: Supression type object; if nothing, a default supression type is constructed
- `m::Int=10`: Krylov subspace dimension, 15 seems to work well for standard grids. Larger m gives slower sampling but better numerical precision.
- `ϵ::Float64=1e-6`: Tykhonoff regularisation
- `n_samples::Int=1000`: Number of MCMC samples to generate
- `burnin::Int=200`: Number of initial samples to discard as burnin
- `thinning::Int=50`: Interval for thinning samples to reduce autocorrelation

# Returns
- A tuple containing:
    - `analysis`: DataFrame with FUBAR analysis results
    - `θ`: Thinned transformed grid samples from the chain

# Description

This function implements a Smooth Kernel Bayesian Density Inference approach to FUBAR analysis. 
It defines a Gaussian model based on the grid, samples from this model using MCMC, 
and processes the samples to generate posterior probabilities of selection.

If no supression type is provided, a default one is constructed based on the grid dimensions
with a fifth degree polynomial is used"""
function FUBAR_analysis(method::SKBDIFUBAR, grid::FUBARGrid{T}; 
    analysis_name = "skbdi_fubar_analysis", 
    volume_scaling = 1.0,
    write = true,
    verbosity = 1,
    posterior_threshold = 0.95,
    distance_function = standard_fubar_distance_function,
    kernel_function = (d, c) -> exp(-d / c^2),
    kernel_parameter_dimension = 1,
    supression_type = nothing,
    m = 10, 
    ϵ = 1e-6, 
    n_samples = 1000, 
    burnin = 200,
    thinning = 50) where {T}

    @assert n_samples > burnin
    if isnothing(supression_type)
        grid_dimension = length(grid.grid_values)
        last_lower_triangular_index = grid_dimension * (grid_dimension - 1) / 2
        supression_dimensions = accumulate((x, i) -> x + (grid_dimension - i),
            0:(grid_dimension-1);
            init=last_lower_triangular_index)
        alphas = 0.1 .* [i for i in 0:(length(supression_dimensions)-1)]
        betas = 0.1 .* [i for i in 1:(length(supression_dimensions))]

        supression_type = SupressionType{T}(alphas, betas,
            supression_dimensions,
            quintic_smooth_transition)
    end
    
    model = define_gaussian_model(grid, 
        distance_function = distance_function,
        kernel_function = kernel_function,
        kernel_parameter_dimension = kernel_parameter_dimension,
        supression_type = supression_type)

    samples, kernel_samples = sample_gaussian_model(model, 
        m = m, 
        n_samples = n_samples,
        burnin = burnin)

    θ = gaussian_sample_postprocessing(model, samples; 
        thinning = thinning, 
        m = m)

    results = FUBAR_bayesian_postprocessing(θ, grid, kernel_samples)
    analysis = FUBAR_tabulate_results(method, results,grid,analysis_name = analysis_name, write = write)
    FUBAR_plot_results(PlotsExtDummy(), method, results, grid, analysis_name = analysis_name, write = write)   
    return analysis, (θ = θ, )
end

function FUBAR_tabulate_results(method::SKBDIFUBAR,results::BayesianFUBARResults, grid::FUBARGrid; analysis_name = "skbdi_fubar_analysis", write = true)
    return FUBAR_tabulate_results(DefaultBayesianFUBARMethod(), results,grid, analysis_name = analysis_name, write = write)
end

## HERE ENDS INTEGRATION WITH THE FUBAR INTERACE

function kernel_matrix_c(grid::FUBARgrid, c::Real)
    n_points = length(grid.alpha_ind_vec)
    K = zeros(eltype(grid.grid_values), n_points, n_points)
    inv_2c_squared = 1.0 / (2 * c^2)  # Precompute this value

    @inbounds for i in 1:n_points
        for j in 1:n_points
            # Calculate distance using alpha_ind_vec and beta_ind_vec like in kernel_matrix
            distance = (grid.alpha_ind_vec[i] - grid.alpha_ind_vec[j])^2 +
                       (grid.beta_ind_vec[i] - grid.beta_ind_vec[j])^2
            K[i, j] = exp(-distance * inv_2c_squared)
        end
    end
    return K
end

function rearrange_kernel_matrix(Σ)
    fubar_to_ess = get_fubar_to_ess_permutation(Int64(sqrt(size(Σ)[1])))
    return Σ[fubar_to_ess, fubar_to_ess]
end

function get_fubar_to_ess_permutation(N::Int)
    ess = generate_ess_indices(N)
    fubar = generate_fubar_indices(N)

    # Create mapping from FUBAR indices to ESS indices
    n_elements = N * N
    perm = zeros(Int, n_elements)

    for i in 1:N
        for j in 1:N
            ess_idx = ess[i, j]
            fubar_idx = fubar[i, j]
            # We want to map FROM fubar TO ess, so:
            perm[ess_idx] = fubar_idx  # This line was wrong before!
        end
    end

    return perm
end

function get_cholesky_total_kernel_matrix(model, gridsize, c, supression_σ; ϵ=1e-6)
    Σ_grid = kernel_matrix_c(model.grid, c)
    Σ_total = zeros(gridsize + 1, gridsize + 1)
    Σ_total[1:gridsize, 1:gridsize] = rearrange_kernel_matrix(Σ_grid)
    Σ_total[gridsize + 1, gridsize + 1] = supression_σ
    return cholesky(Symmetric(Σ_total + ϵ * I)).L
end

struct MatrixInterpolatingGPFUBAR
    model::RJGPModel
    αs::Vector{Float64}
    βs::Vector{Float64}
    dimension::Int64
    fubar_to_ess_perm::Vector{Int64}
    pre_determined_c::Vector{Float64}
    pre_computed_cholesky_matrices::Vector{Matrix{Float64}}
    kernel_matrix_indicator_prior::Distribution
    supression_σ::Float64
    function MatrixInterpolatingGPFUBAR(model::RJGPModel, αs::Vector{Float64}, βs::Vector{Float64}, dimension::Int64, fubar_to_ess_perm::Vector{Int64}, pre_determined_c::Vector{Float64}, supression_σ::Float64)
        pre_computed_cholesky_matrices = [get_cholesky_total_kernel_matrix(model, dimension, c, supression_σ) for c in pre_determined_c]
        kernel_matrix_indicator_prior = Normal(Float64(length(pre_determined_c)/2), Float64(length(pre_determined_c)/4)) # Interesting
        return new(model, αs, βs, dimension, fubar_to_ess_perm, pre_determined_c, pre_computed_cholesky_matrices, kernel_matrix_indicator_prior, supression_σ)
    end
end

function interpolate_matricies(A, B, t)
    return A * (1 - t) + B * t
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


function supress_vector(θ, αs, βs, dimensions)
    T = eltype(θ)
    grid_θ = @view θ[1:end-1]  # Use view instead of copying
    y = θ[end]

    s = softmax(grid_θ)
    M = ones(T, length(s))

    @inbounds for i in 2:length(dimensions)
        start_index = dimensions[i-1] + 1
        end_index = dimensions[i]
        f_i = quintic_smooth_transition(y, αs[i], βs[i])
        M[start_index:end_index] .= f_i
    end

    A = s .* M
    T_sum = sum(A)
    return A ./ T_sum
end

function supression_loglikelihood(model::RJGPModel, θ, αs, βs, dimensions)
    T = eltype(θ)
    softmax_grid_θ = supress_vector(θ, αs, βs, dimensions)
    return sum(log.(softmax_grid_θ[model.ess_to_fubar_perm]' * model.grid.cond_lik_matrix))
end

function get_interpolated_kernel_matrix(problem::MatrixInterpolatingGPFUBAR, kernel_matrix_indicator)
    lower_matrix = Int64(max(1,floor(kernel_matrix_indicator)))
    upper_matrix = Int64(max(2,min(length(problem.pre_computed_cholesky_matrices),ceil(kernel_matrix_indicator))))
    t = kernel_matrix_indicator - lower_matrix
    return interpolate_matricies(problem.pre_computed_cholesky_matrices[lower_matrix], problem.pre_computed_cholesky_matrices[upper_matrix], t)
end

function mvnormal_logpdf_cholesky(x, μ, L)
    n = length(μ)
    
    # Compute log determinant of Σ using L
    logdet_Σ = 2 * sum(log.(diag(L)))
    
    # Compute (x-μ)ᵀ × Σ^(-1) × (x-μ) using Cholesky
    diff = x - μ
    z = L \ diff  # Solves L × z = diff
    quad = dot(z, z)
    
    # Compute log PDF
    logpdf = -0.5 * (n * log(2π) + logdet_Σ + quad)

    return logpdf
end

function (problem::MatrixInterpolatingGPFUBAR)(θ)
    kernel_matrix_indicator = θ[end] # This is the parameter that controls which kernel matrix to use
    kernel_matrix = get_interpolated_kernel_matrix(problem, kernel_matrix_indicator)
    log_likelihood = supression_loglikelihood(problem.model, θ[1:(end-1)], problem.αs, problem.βs, problem.dimension)
    kernel_matrix_indicator_log_prior = logpdf(problem.kernel_matrix_indicator_prior, kernel_matrix_indicator)    
    grid_log_prior = mvnormal_logpdf_cholesky(θ[1:(end-1)], zeros(length(θ[1:(end-1)])), kernel_matrix)
    return log_likelihood + kernel_matrix_indicator_log_prior + grid_log_prior
end

LogDensityProblems.logdensity(problem::MatrixInterpolatingGPFUBAR, θ) = problem(θ)
LogDensityProblems.dimension(problem::MatrixInterpolatingGPFUBAR) = problem.model.dimension + 2
LogDensityProblems.capabilities(::Type{MatrixInterpolatingGPFUBAR}) = LogDensityProblems.LogDensityOrder{0}()

function matrix_interpolating_gp_fubar_HMC_sample(model::RJGPModel, n_samples::Int64)
    dimensions = accumulate((x, i) -> x + (20 - i), 0:(20-1); init=190) #Hard coded 20,190 for now
    αs = 0.1 .* [i for i in 0:(length(dimensions)-1)]
    βs = 0.1 .* [i for i in 1:(length(dimensions))]
    # Pre-compute the permutation
    fubar_to_ess_perm = get_fubar_to_ess_permutation(Int64(sqrt(model.dimension)))

    problem = MatrixInterpolatingGPFUBAR(model, αs, βs, 400, fubar_to_ess_perm, collect(0.5:0.01:2.0), 1.0)


    model = AdvancedHMC.LogDensityModel(LogDensityProblemsAD.ADgradient(Val(:Zygote), problem))
    δ = 0.8
    sampler = NUTS(δ)
    samples = AbstractMCMC.sample(
        model,
        sampler,
        n_samples;
        initial_params=randn(problem.model.dimension + 2),
    )
    return samples
end
using  Statistics, Distributions, MCMCChains, LinearAlgebra, Phylo, FASTX, MolecularEvolution, CodonMolecularEvolution, Plots, EllipticalSliceSampling, AbstractMCMC, NNlib, StatsBase

#= This struct is a bit f****d
    because it is only applicable to difFUBAR
    even though the suppression scheme is a lot more general...
    I plan on fixing this later.
    
    Maybe I should try to divide it into
    unspecified suppression parameters
    unspecified unsuppressed parameters
    and then have another struct for the ESS model
    which adds additional parameters like
    sigma_c or Sigma_c stdev/covmat where c are the "kernel parameters"
    sigma_s or Sigma_s stdev/covmat where s are the suppression parameters
    and the kernel function 
        
    I WILL do this as soon as I have a working version of the code =#

struct SuppressionGrid
    # Array of precomputed likelihood values 
    # indexed by category, site
    con_lik_matrix::Array{Float64}
    # Array of precomputed log-likelihood values 
    # indexed by category, site
    log_con_lik_matrix::Array{Float64}
    # Indexed codon, param where 
    # 1, 2, 3, 4 mean alpha, omega_1, omega_2 and omega_BG, respectvely
    codon_param_mat::Array{Float64}
    codon_param_index_mat::Array{Int} # Indexed codon, param index in the grids
    alpha_grid::Vector{Float64} # Vector of values of alpha 
    omega_grid::Vector{Float64} # Vector of values of omega
    background_omega_grid::Vector{Float64} # Vector of values of background omega
    param_kinds::Vector{String} # Vector of parameter names
    suppression_kinds::Vector{String} # Vector of suppression parameter names
    transition_functions::Vector{Function} # Vector of transition functions for the suppression parameters
    # Array of boolean masks for suppression for parameter values 
    # indexed suppression_parameter, parameter, grid_index 
    # where the order is to be the same as in param_kinds, suppression_kinds
    masks::Array{Bool}
end

# It makes sense to reparametrize to logit(u_theta) the entire way...
function probability_vector_given_parameters(grid::SuppressionGrid, parameters::Vector{Float64}, sigma_s::Float64, sigma_c::Float64, kernel_function::Function, square_distance_matrix::Matrix{Int64}; epsilon::Float64=1e-6)
    # Parameters divided into 
    # 1. kernel parameters, 2. suppression parameters, 3. unsuppressed parameters
    # Now only one suppression parameter c is used...
    # sqrt(Sigma) * parameters = vcat((c, log s, "unsoftmaxed" unsuppressed parameters))
    num_suppression_parameters = size(grid.masks)[1]
    num_kernel_parameters = length(parameters) - num_suppression_parameters - size(grid.log_con_lik_matrix)[1]
    suppression_parameters = sigma_s .* parameters[num_kernel_parameters + 1:num_suppression_parameters + num_kernel_parameters]
    unsuppressed_cov_mat = kernel_function(sigma_c * parameters[1:num_kernel_parameters], square_distance_matrix)
    probability_vector = softmax(CodonMolecularEvolution.krylov_sqrt_times_vector(unsuppressed_cov_mat + epsilon * I, parameters[num_kernel_parameters + num_suppression_parameters+1:end]))
    for i = 1:num_suppression_parameters
        probability_vector[grid.masks[i, :]] .*= grid.transition_functions[i](suppression_parameters[i])
    end
    return probability_vector / sum(probability_vector) # Normalization
end

function log_likelihood(grid::SuppressionGrid, parameters::Vector{Float64}, sigma_c::Float64, sigma_s::Float64, kernel_function::Function, square_distance_matrix::Matrix{Int64}; epsilon::Float64=1e-6)
    # l(D|unsuppressed_theta, s) = ∑_i log(∑_k P(site i|C_k) P(C_k|unsuppressed_theta, s))
    return sum(log.(grid.con_lik_matrix' * probability_vector_given_parameters(grid, parameters, sigma_s, sigma_c, kernel_function, square_distance_matrix, epsilon=epsilon)))
end

function define_ess_model(grid::SuppressionGrid, num_kernel_parameters::Int64, kernel_function::Function, sigma_s::Float64, sigma_c::Float64, square_distance_matrix::Matrix{Int64})
    """
    # Arguments
    - `grid::FUBARGrid{T}`: Grid to perform inference on
    - `sigma_s::Float64`:Standard deviation for the normal prior for suppression parameters
    - `Sigma_0::Float64`:Covariance matrix for the mvnormal prior for logit unsuppressed parameters
    TODO: Add the parameter c for the kernel function as "learnable"
    """
    num_s = size(grid.masks)[1]
    num_u_theta = size(grid.log_con_lik_matrix)[1]
    total_dimension = num_kernel_parameters + num_s + num_u_theta
    #= s_cov_mat = diagm(sigma_s^2 .* ones(Float64, num_s))
    prior_cov_mat = [s_cov_mat zeros(Float64, (num_s, num_u_theta));
        zeros(Float64, (num_u_theta, num_s)) Sigma_0] =#
    prior = MvNormal(zeros(Float64, total_dimension), I)
    log_likelihood_one_arg = parameters -> log_likelihood(grid, parameters, sigma_c, sigma_s, kernel_function, square_distance_matrix)
    return ESSModel(prior, log_likelihood_one_arg)
end

function calculate_alloc_grid_and_theta(grid::SuppressionGrid, ambient_samples::Vector{Vector{Float64}}, burnin::Int=0)
    """
    Samples the allocation grid from the ambient samples
    TODO: WHY DOES THIS TAKE SO MUCH TIME x_x
    """
    n_samples = size(ambient_samples, 1)
    n_categories = size(grid.con_lik_matrix, 1)
    n_sites = size(grid.con_lik_matrix, 2)
    alloc_grid = zeros(Int64, size(grid.con_lik_matrix))
    theta = zeros(Float64, n_categories)
    for i = burnin+1:size(ambient_samples, 1)
        probability_vector = probability_vector_given_parameters(grid, ambient_samples[i])
        theta .+= probability_vector # Accumulate the probability vector
        v = zeros(Float64, n_categories)
        for i = 1:n_sites
            for j = 1:n_categories
                v[j] = probability_vector[j] * grid.con_lik_matrix[j, i]
            end
            samp = sample(1:n_categories, Weights(v))
            alloc_grid[samp, i] += 1 # Increment the allocation for the sampled category
        end
    end
    return alloc_grid, theta ./ (n_samples - burnin)
end

function fast_cov_mat_hedwigs_kernel(c::Vector{Float64}, square_distance_matrix::Matrix{Int64})
    return exp(-1 / c[1]^2).^square_distance_matrix
end

function generate_square_l2_distance_matrix(codon_param_index_mat::Array{Int})
    # Generate a distance matrix for the codon parameters
    n = size(codon_param_index_mat, 1)
    square_distance_matrix = zeros(Int64, n, n)
    for i = 1:n
        for j = 1:i
            square_distance_matrix[i, j] = sum((codon_param_index_mat[i, :] .- codon_param_index_mat[j, :]).^2)
            square_distance_matrix[j, i] = square_distance_matrix[i, j] # Symmetric matrix
        end
    end
    return square_distance_matrix
end

function skbdifFUBAR(seqnames, seqs, treestring, tags, outpath,
    tag_colors; pos_thresh=0.95, iters=2500,
    burnin::Int=div(iters, 5), concentration=0.1, binarize=false, verbosity=1,
    exports=true, exports2json=false, code=MolecularEvolution.universal_code,
    optimize_branch_lengths=false, version=nothing, t=0)
    # TODO: Split this up into multiple fcts because the setup is mostly shared with difFUBAR
    total_time = @elapsed begin
    analysis_name = outpath
    leaf_name_transform = CodonMolecularEvolution.generate_tag_stripper(tags)
    plot_collection = NamedTuple[]
    tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name,
        treestring,
        tags,
        tag_colors=tag_colors,
        exports=exports,
        verbosity=verbosity,
        disable_binarize=!binarize, plot_collection=plot_collection)
    ((tree, LL, alpha, beta, GTRmat, F3x4_freqs, eq_freqs), fit_time) =
        @timed CodonMolecularEvolution.difFUBAR_global_fit_2steps(seqnames,
            seqs,
            tree,
            leaf_name_transform,
            code,
            verbosity=verbosity,
            optimize_branch_lengths=optimize_branch_lengths)

    ((con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds, shallow_tree, background_omega_grid, codon_param_index_vec), grid_time) =
        @timed CodonMolecularEvolution.difFUBAR_grid(tree,
            tags,
            GTRmat,
            F3x4_freqs,
            code,
            verbosity=verbosity,
            foreground_grid=4,
            background_grid=2,
            version=version,
            t=t)

    transition_functions = [s -> CodonMolecularEvolution.quintic_smooth_transition(s, 0, 1) for i = 1:4]
    # Define the masks for the suppression parameters
    masks = ones(Bool, (4, size(con_lik_matrix)[1]))
    
    codon_param_mat = reduce(hcat, codon_param_vec)'#reshape(codon_param_vec, (size(codon_param_vec, 1), 3)) # Reshape to codon, param
    codon_param_index_mat = reduce(hcat, codon_param_index_vec)' # Reshape to codon, param index

    masks[1, :] .= codon_param_mat[:,2] .> 1 # ω_1 > 1
    masks[2, :] .= codon_param_mat[:,3] .> 1 # ω_2 > 1
    masks[3, :] .= codon_param_mat[:,2] .> codon_param_mat[:,3] # ω_1 > ω_2
    masks[4, :] .= codon_param_mat[:,3] .> codon_param_mat[:,2] # ω_1 < ω_2
    grid = SuppressionGrid(con_lik_matrix, log_con_lik_matrix, codon_param_mat, codon_param_index_mat, alphagrid, omegagrid, background_omega_grid, param_kinds, ["ω_1>ω_2", "ω_1<ω_2", "ω_1>1", "ω_2>1"], transition_functions, masks)
    #squared_distance_function = (i, j) -> sum((grid.codon_param_index_mat[i, :].- grid.codon_param_index_mat[j, :]).^2) # Squared distance function for the kernel
    #c = 1.0
    square_distance_matrix = generate_square_l2_distance_matrix(grid.codon_param_index_mat)
    #kernel_function = (i, j) -> exp(-squared_distance_function(i, j) / c) # Gaussian kernel
    #sigma_0 = [kernel_function(i, j) for i in 1:size(codon_param_mat, 1), j in 1:size(codon_param_mat, 1)]
    sigma_s = 0.1 # Standard deviation for the normal prior for suppression parameters
    sigma_c = 1.0 # Standard deviation for the normal prior for kernel parameter
    square_distance_matrix = generate_square_l2_distance_matrix(grid.codon_param_index_mat)
    ess_model = define_ess_model(grid, 1, fast_cov_mat_hedwigs_kernel, sigma_s, sigma_c, square_distance_matrix) # Define the ESS model
    # Run the ESS sampling
    sample_time = @elapsed begin
    ambient_samples = AbstractMCMC.sample(ess_model, ESS(), iters,
        progress=true)
    alloc_grid, theta = calculate_alloc_grid_and_theta(grid, ambient_samples, burnin)
    end # End of sample_time
    df, plots_named_tuple = CodonMolecularEvolution.difFUBAR_tabulate_and_plot(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors, verbosity=verbosity, exports=exports)
    end # End of total_time
    # TODO: Ask Maxi why I don't get all of the plots that difFUBAR does ????????? :c
    json = CodonMolecularEvolution.dNdS2JSON(CodonMolecularEvolution.difFUBAR2JSON(), (outpath=analysis_name, df=df, θ=theta, posterior_mat=alloc_grid ./ sum(alloc_grid, dims=1), categories=reduce(hcat, codon_param_vec)', tags=tags, tree=shallow_tree, LL=LL, timers = (total_time, fit_time, grid_time, sample_time), treestring=treestring, seqnames=seqnames, seqs=seqs, leaf_name_transform=leaf_name_transform, pos_thresh=pos_thresh, iters=iters, burnin=burnin, concentration=concentration, binarize=binarize, exports=exports2json))
    push!(plot_collection, plots_named_tuple)
    return df, (alloc_grid, ambient_samples, grid, tag_colors), merge(plot_collection...)
end


analysis_name = "output/Ace2"
seqnames,seqs = read_fasta("test/data/Ace2_tiny/Ace2_tiny_tagged.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("test/data/Ace2_tiny/Ace2_tiny_tagged.nex")
df, results = skbdifFUBAR(seqnames, seqs, treestring, tags, analysis_name, tag_colors)








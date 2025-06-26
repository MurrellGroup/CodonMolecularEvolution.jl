using  Distributions, MCMCChains, LinearAlgebra, Phylo, FASTX, MolecularEvolution, CodonMolecularEvolution, Phylo, EllipticalSliceSampling

#= This struct is a bit f****d
    because it is only applicable to difFUBAR
    even though the suppression scheme is a lot more general...
    I plan on fixing this later.=#
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
    alpha_grid::Vector{Float64} # Vector of values of alpha 
    omega_grid::Vector{Float64} # Vector of values of omega
    param_kinds::Vector{String} # Vector of parameter names
    suppression_kinds::Vector{String} # Vector of suppression parameter names
    transition_functions::Vector{Function} # Vector of transition functions for the suppression parameters
    # Array of boolean masks for suppression for parameter values 
    # indexed suppression_parameter, parameter, grid_index 
    # where the order is to be the same as in param_kinds, suppression_kinds
    masks::Array{Bool}
end

# It makes sense to reparametrize to logit(u_theta) the entire way...
function probability_vector_given_parameters(grid::SuppressionGrid, parameters::Vector{Float64})
    num_suppression_parameters = size(grid.masks)[1]
    suppression_parameters = parameters[1:num_suppression_parameters]
    probability_vector = softmax(parameters[num_suppression_parameters+1:end]) # TODO: Find out if using softmax is correct? Feels like should be expit but idkkkk
    for i = 1:num_suppression_parameters
        probability_vector[grid.masks[i]] .*= grid.transition_functions[i](suppression_parameters[i])
    end
    return probability_vector / sum(probability_vector) # Normalization
end

function log_likelihood(grid::SuppressionGrid, parameters::Vector{Float64})
    # l(D|unsuppressed_theta, s) = ∑_i log(∑_k P(site i|C_k) P(C_k|unsuppressed_theta, s))
    return sum(log.(grid.con_lik_matrix' * probability_vector_given_parameters(grid, parameters)))
end



function define_ess_model(grid::SuppressionGrid, sigma_s::Float64, Sigma_0::Matrix{Float64})
    """
    # Arguments
    - `grid::FUBARGrid{T}`: Grid to perform inference on
    - `sigma_s::Float64`:Standard deviation for the normal prior for suppression parameters
    - `Sigma_0::Float64`:Covariance matrix for the mvnormal prior for logit unsuppressed parameters
    """
    num_s = size(grid.masks)[1]
    num_u_theta = size(grid.log_con_lik_matrix)[1]
    s_cov_mat = diagm(sigma_s^2 .* ones(Float64, num_s))
    prior_cov_mat = [s_cov_mat zeros(Float64, (num_s, num_u_theta));
        zeros(Float64, (num_u_theta, num_s)) Sigma_0]
    prior = MvNormal(zeros(Float64, num_s + num_u_theta), prior_cov_mat)
    log_likelihood = parameters -> log_likelihood(grid, parameters)
    return ESSModel(prior, log_likelihood)
end

function distance_grid_index(grid::SuppressionGrid, i::Int, j::Int)
    """
    Returns the distance between two grid indices i and j
    """
    index_space_vectors = zeros(Float64, (2, 4))
    indices = [i, j]
    for k = 1:2
        index_space_vectors[k, 1] = findfirst(==(grid.codon_param_mat[indices[k], 1]), grid.alpha_grid) # Find the index of alpha
        for l = 2:4
            index_space_vectors[k, l] = findfirst(==(grid.codon_param_mat[indices[k], l]), grid.omega_grid) # Find the index of omega
        end
    end
    return sum((index_space_vectors[1, :] .- index_space_vectors[2, :]).^2) # Squared Euclidean distance in index space
end

function skbdifFUBAR(seqnames, seqs, treestring, tags, outpath,
    tag_colors; pos_thresh=0.95, iters=2500,
    burnin::Int=div(iters, 5), concentration=0.1, binarize=false, verbosity=1,
    exports=true, exports2json=false, code=MolecularEvolution.universal_code,
    optimize_branch_lengths=false, version=nothing, t=0)
    # TODO: Split this up into multiple fcts because the setup is mostly shared with difFUBAR

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

    ((con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _, shallow_tree), grid_time) =
        @timed CodonMolecularEvolution.difFUBAR_grid(tree,
            tags,
            GTRmat,
            F3x4_freqs,
            code,
            verbosity=verbosity,
            foreground_grid=2,
            background_grid=2,
            version=version,
            t=t)

    log_con_lik_matrix = log.(con_lik_matrix)
    transition_functions = [s -> CodonMolecularEvolution.quintic_smooth_transition(s, 0, 1) for i = 1:4]
    # Define the masks for the suppression parameters
    masks = ones(Bool, (4, size(con_lik_matrix)[1]))
    codon_param_mat = reduce(hcat, codon_param_vec)'#reshape(codon_param_vec, (size(codon_param_vec, 1), 3)) # Reshape to codon, param


    masks[1, :] .= codon_param_mat[:,2] .> 1 # ω_1 > 1
    masks[2, :] .= codon_param_mat[:,3] .> 1 # ω_2 > 1
    masks[3, :] .= codon_param_mat[:,2] .> codon_param_mat[:,3] # ω_1 > ω_2
    masks[4, :] .= codon_param_mat[:,3] .> codon_param_mat[:,2] # ω_1 < ω_2
    grid = SuppressionGrid(con_lik_matrix, log_con_lik_matrix, codon_param_mat, alphagrid, omegagrid, ["α,ω_1,ω_2,ω_BG"], ["ω_1>ω_2", "ω_1<ω_2", "ω_1>1", "ω_2>1"], transition_functions, masks)

    squared_distance_function = (i, j) -> distance_grid_index(grid, i, j)#sum((i .- j).^2) # Squared Euclidean distance
    c = 1.0
    kernel_function = (i, j) -> exp(-squared_distance_function(i, j) / c) # Gaussian kernel
    sigma_0 = [kernel_function(i, j) for i in 1:size(codon_param_mat, 1), j in 1:size(codon_param_mat, 1)]
    sigma_s = 0.1 # Standard deviation for the normal prior for suppression parameters
    ess_model = define_ess_model(grid, sigma_s, sigma_0)

    # TODO: Sample from the posterior
    # TODO: summary statistics, visualizations, save results
end


analysis_name = "output/Ace2"
seqnames,seqs = read_fasta("test/data/Ace2_tiny/Ace2_tiny_tagged.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("test/data/Ace2_tiny/Ace2_tiny_tagged.nex")
skbdifFUBAR(seqnames, seqs, treestring, tags, analysis_name, tag_colors)








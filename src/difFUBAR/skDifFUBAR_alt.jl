using  Statistics, Distributions, MCMCChains, LinearAlgebra, Phylo, FASTX, MolecularEvolution, CodonMolecularEvolution, Plots, EllipticalSliceSampling, AbstractMCMC, NNlib, StatsBase

"""
SKBDIModel: holds the model parameters for
a general skbdi model.
parameter_grids: A vector of vectors, where each inner vector contains 
                 the parameter values for a specific grid.
parameter_names: A vector of strings, where each string is the name 
                 of a parameter corresponding to the parameter values in parameter_grids.
masks: A matrix of booleans with suppress[suppresion_parameter, category] = true 
       if the category is suppressed by the suppression parameter.
transition_functions: A vector of transition functions F:R -> [0, 1] for the skbdi model
log_con_lik_matrix: A matrix of log likelihoods for each category.
con_lik_matrix: A matrix of likelihoods for each category.
codon_param_vec: A vector of vectors indexed category, parameter containing
                 parameter values for each category
codon_param_index_vec: A vector of integers indexed category, parameter containing the indices 
                       in the parameter grids for each category.
ambient_to_parameter_transform: A function that transforms an ambient sample into the parameter space.
kernel_dim: The number of kernel parameters.
suppression_dim: The number of suppression parameters.
unsuppressed_dim: The number of unsuppressed parameters.
total_dim: The total number of parameters in the model.

# TODO: Only ESS for now, make it more general later.
"""
struct SKBDIModel
    # Model Parameters:
    parameter_grids::Vector{Vector{Float64}}
    parameter_names::Vector{String}
    masks::Matrix{Bool}
    transition_functions::Vector{Function}
    log_con_lik_matrix::Matrix{Float64}
    con_lik_matrix::Matrix{Float64}
    codon_param_vec::Vector{Vector{Float64}}
    codon_param_index_vec::Vector{Int64}
    ambient_to_parameter_transform::Function{Vector{Float64}, Vector{Float64}}
    kernel_dim::Int64
    suppression_dim = size(masks)[1]
    unsuppressed_dim = size(log_con_lik_matrix)[1]
    total_dim = kernel_dim + suppression_dim + unsuppressed_dim
    # MCMC Parameters maybe move these to a separate struct?
    #= sampler::Function
    online::Bool =#
    num_samples::Int64
    burnin::Int64
end

function split_parameters(model::SKBDIModel, parameters::Vector{Float64})
    """
    Splits the parameters into kernel, suppression, and unsuppressed parameters.
    """
    kernel_parameters = parameters[1:model.kernel_dim]
    suppression_parameters = parameters[model.kernel_dim + 1:model.kernel_dim + model.suppression_dim]
    unsuppressed_parameters = parameters[model.kernel_dim + model.suppression_dim + 1:end]
    return kernel_parameters, suppression_parameters, unsuppressed_parameters
end

function calculate_probability_vector(model::SKBDIModel, parameters::Vector{Float64})
    """
    Computes the probability vector for the model given the parameters.
    parameters start with kernel parameters, followed by suppression parameters, 
    and then unsuppressed parameters.
    """
    _, suppression_parameters, unsuppressed_parameters = split_parameters(model, parameters)
    probability_vector = softmax(unsuppressed_parameters)
    # Apply the transition functions for each suppression parameter
    for i in 1:model.suppression_dim
        probability_vector[model.masks[i, :]] .*= model.transition_functions[i](suppression_parameters)
    end
    return probability_vector ./ sum(probability_vector)
end

function log_likelihood(model::SKBDIModel, ambient_sample::Vector{Float64})
    """
    Computes the log-likelihood of the model given the parameters.
    """
    parameters = model.ambient_to_parameter_transform(ambient_sample)
    probability_vector = calculate_probability_vector(model, parameters)
    return sum(log.(model.con_lik_matrix' * probability_vector))
end

function define_ess_model(model::SKBDIModel)
    """
    a wrapper function to define the model for ESS sampling.
    It takes a SKBDIModel and returns an EllipticalSliceSampling.Model.
    """
    return ESSModel(
        MvNormal(zeros(Float64, model.total_dim), I),
        (ambient_sample) -> log_likelihood(model, ambient_sample),
    )
end

function hedwigs_ambient_to_parameter_transform(model::SKBDIModel, ambient_sample::Vector{Float64}, kernel_stddev::Float64, suppression_stddev::Float64, square_distance_matrix::Matrix{Int64}, kernel_function::Function, epsilon::Float64=1e-6)
    """
    Transforms an ambient sample (~N(0, I)) into the parameter space (~N(0, Sigma)).
    """
    kernel_parameters, suppression_parameters, unsuppressed_parameters = split_parameters(model, ambient_sample)
    kernel_parameters = kernel_stddev * kernel_parameters
    suppression_parameters = suppression_stddev * suppression_parameters
    covariance_matrix = kernel_function(kernel_parameters, square_distance_matrix) + epsilon * I # Tykhonoff regularization
    unsuppressed_parameters = CodonMolecularEvolution.krylov_sqrt_times_vector(covariance_matrix, 
        kernel_parameters)
    return vcat(kernel_parameters, suppression_parameters, unsuppressed_parameters)
end

function fast_cov_mat_hedwigs_kernel(c::Vector{Float64}, square_distance_matrix::Matrix{Int64})
    """
    A fast implementation of Hedwig's kernel for the covariance matrix.
    """
    return exp(-c[1]^2).^square_distance_matrix
end

function generate_square_l2_distance_matrix(codon_param_index_vec::Vector{Int64})
    """
    Generates a square matrix of L2 distances between the categories.
    """
    num_categories = length(codon_param_index_vec)
    square_distance_matrix = zeros(Int64, num_categories, num_categories)
    for i in 1:num_categories
        for j in 1:i
            # Note: The following fails if i, j are not the same length.
            square_distance_matrix[i, j] = sum((codon_param_index_vec[i] .- codon_param_index_vec[j]).^2)
            square_distance_matrix[j, i] = square_distance_matrix[i, j]  # Symmetric matrix
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
    masks[1, :] = [c[2] > 1 for c in codon_param_vec] # omega_1 > 1
    masks[2, :] = [c[3] > 1 for c in codon_param_vec] # omega_2 > 1
    masks[3, :] = [c[2] > c[3] for c in codon_param_vec] # omega_1 > omega_2
    masks[4, :] = [c[3] > c[2] for c in codon_param_vec] # omega_2 > omega_1
    
    square_distance_matrix = generate_square_l2_distance_matrix(codon_param_index_vec)
    kernel_stddev = 0.1 # example values idk what these should be xD
    suppression_stddev = 0.1
    model = SKBDIModel(
        [alphagrid, omegagrid, background_omega_grid],
        ["alpha", "omega_1", "omega_2", "background_omega"],
        masks,
        transition_functions,
        log_con_lik_matrix,
        con_lik_matrix,
        codon_param_vec,
        codon_param_index_vec,
        (ambient_sample, c) -> hedwigs_ambient_to_parameter_transform(model, ambient_sample), # ambient_to_parameter_transform will be defined later
        length(alphagrid), # kernel_dim
    )


    
    
    return df, (alloc_grid, ambient_samples, grid, tag_colors), merge(plot_collection...)
end



using Statistics, Combinatorics, Distributions, MCMCChains, LinearAlgebra, Phylo, FASTX, MolecularEvolution, CodonMolecularEvolution, Plots, EllipticalSliceSampling, AbstractMCMC, NNlib, StatsBase, ProgressMeter, JLD2

"""
# GeneralizedFUBARModel
A general framework for the class of models to which FUBAR and skbdi belong.
## Fields
- n_parameters: Number of parameters in the model
- n_categories: Number of categories in the model
- log_likelihood: log likelihood for a probability vector for the categories
- to_probability_vector: returns the probability vector for a given set of parameters
- prior: the prior distribution in parameter space.
- con_lik_matrix: the matrix of conditional distributions for the sites.
"""
struct GeneralizedFUBARModel
    n_parameters::Int64
    n_categories::Int64
    log_likelihood::Function
    to_probability_vector::Function
    prior::Distribution
    con_lik_matrix::Matrix{Float64}
end

function log_posterior(model::GeneralizedFUBARModel, parameters::Vector{Float64})
    return model.log_likelihood(model.parameters_to_probability_vector(parameters)) + model.log_prior(parameters)
end

"""
# SKBDIModel: 
holds the model parameters for
a general skbdi model.
## Fields:
- parameter_grids: A vector of vectors, where each inner vector contains 
                 the parameter values for a specific grid.
- parameter_names: A vector of strings, where each string is the name 
                 of a parameter corresponding to the parameter values in parameter_grids.
- hypothesis_masks: A matrix of booleans with mask[hypothesis, category] = true 
       if the hypothesis is true for the category.
- transition_functions: A vector of transition functions F:R -> [0, 1] for the skbdi model
- log_con_lik_matrix: A matrix of log likelihoods for each category.
- con_lik_matrix: A matrix of likelihoods for each category.
- codon_param_vec: A vector of vectors indexed category, parameter containing
                 parameter values for each category
- codon_param_index_vec: A vector of integers indexed category, parameter containing the indices 
                       in the parameter grids for each category.
- ambient_to_parameter_transform: A function that transforms an ambient sample into the parameter space.
- kernel_dim: The number of kernel parameters.
- masks: The disjoint masks used by the model when applying suppression parameters.
- suppression_dim: The number of suppression parameters.
- unsuppressed_dim: The number of unsuppressed parameters.
- total_dim: The total number of parameters in the model.
"""
struct SKBDIModel
    # Model Parameters:
    parameter_grids::Vector{Vector{Float64}}
    parameter_names::Vector{String}
    hypothesis_masks::Matrix{Bool}
    transition_function::Function
    log_con_lik_matrix::Matrix{Float64}
    con_lik_matrix::Matrix{Float64}
    codon_param_vec::Vector{Vector{Float64}}
    codon_param_index_vec::Vector{Vector{Int64}}
    ambient_to_parameter_transform::Function
    kernel_dim::Int64
    ## Derived fields:
    masks::Matrix{Bool}
    mask_subset_indicators::Matrix{Bool}
    suppression_dim::Int64
    unsuppressed_dim::Int64
    total_dim::Int64
end

function GeneralizedFUBARModel(model::SKBDIModel)
    return GeneralizedFUBARModel(
        model.total_dim,
        model.unsuppressed_dim,
        a -> log_likelihood(model, a),
        a -> to_probability_vector(model, a),
        MvNormal(zeros(model.total_dim), I),
        model.con_lik_matrix)
end

"""
# generate_mutually_exclusive_masks
generates additional suppression masks for each non-empty combination of hypotheses
"""
function gen_disjoint_masks(hypothesis_masks::Matrix{Bool})
    disjoint_masks = Vector{Bool}[]
    mask_subset_indicators = Vector{Bool}[]
    # Go through subsets (grows super fast with number of hypotheses)
    for i = 1:(2^(size(hypothesis_masks, 1))-1)
        b = [(i >> (j - 1) & 1) == 1 for j = 1:size(hypothesis_masks, 1)]
        # u (⋃_{j:b_j == 1} masks[j])∖(⋂_{j:b_j==0} masks[j])
        u = prod(hypothesis_masks[b, :], dims=1) .& .!(sum(hypothesis_masks[.!b, :], dims=1) .> 0)
        if any(u)
            push!(disjoint_masks, Matrix(u)[1, :]) # Ugh annoying typecasting bitmatrix → vector{bool}
            push!(mask_subset_indicators, b)
        end
    end
    return collect(hcat(disjoint_masks...)'), collect(hcat(mask_subset_indicators...)')
end

function SKBDIModel(parameter_grids::Vector{Vector{Float64}},
    parameter_names::Vector{String},
    hypothesis_masks::Matrix{Bool},
    transition_function::Function,
    log_con_lik_matrix::Matrix{Float64},
    con_lik_matrix::Matrix{Float64},
    codon_param_vec::Vector{Vector{Float64}},
    codon_param_index_vec::Vector{Vector{Int64}},
    ambient_to_parameter_transform::Function,
    kernel_dim::Int64)
    masks, mask_subset_indicators = gen_disjoint_masks(hypothesis_masks)
    suppression_dim = size(masks)[1]
    unsuppressed_dim = size(log_con_lik_matrix)[1]
    total_dim = kernel_dim + suppression_dim + unsuppressed_dim


    return SKBDIModel(parameter_grids,
        parameter_names,
        hypothesis_masks,
        transition_function,
        log_con_lik_matrix,
        con_lik_matrix,
        codon_param_vec,
        codon_param_index_vec,
        ambient_to_parameter_transform,
        kernel_dim,
        masks,
        mask_subset_indicators,
        suppression_dim,
        unsuppressed_dim,
        total_dim)
end

function split_parameters(model::SKBDIModel, parameters::Vector{Float64})
    """
    Splits the parameters into kernel, suppression, and unsuppressed parameters.
    """
    kernel_parameters = parameters[1:model.kernel_dim]
    suppression_parameters = parameters[model.kernel_dim+1:model.kernel_dim+model.suppression_dim]
    unsuppressed_parameters = parameters[model.kernel_dim+model.suppression_dim+1:end]
    return kernel_parameters, suppression_parameters, unsuppressed_parameters
end

function to_probability_vector(model::SKBDIModel, ambient_sample::Vector{Float64})
    """
    Computes the probability vector for the model given the parameters.
    parameters start with kernel parameters, followed by suppression parameters, 
    and then unsuppressed parameters.
    """
    parameters = model.ambient_to_parameter_transform(ambient_sample)
    _, suppression_parameters, unsuppressed_parameters = split_parameters(model, parameters)
    probability_vector = softmax(unsuppressed_parameters)
    # Apply the transition functions for each suppression parameter
    for i in 1:model.suppression_dim
        probability_vector[model.masks[i, :]] .*= model.transition_function(suppression_parameters[i])
    end
    return probability_vector ./ sum(probability_vector)
end

function log_likelihood(model::SKBDIModel, ambient_sample::Vector{Float64})
    """
    Computes the log-likelihood of the model given the parameters.
    """
    probability_vector = to_probability_vector(model, ambient_sample)
    return sum(log.(model.con_lik_matrix' * probability_vector))

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
        unsuppressed_parameters)
    return vcat(kernel_parameters, suppression_parameters, unsuppressed_parameters)
end

function fast_cov_mat_hedwigs_kernel(c::Vector{Float64}, square_distance_matrix::Matrix{Int64})
    """
    A fast implementation of Hedwig's kernel for the covariance matrix.
    """
    return exp(-c[1]^2) .^ square_distance_matrix
end

function generate_square_l2_distance_matrix(codon_param_index_vec::Vector{Vector{Int64}})
    """
    Generates a square matrix of L2 distances between the categories.
    """
    num_categories = length(codon_param_index_vec)
    square_distance_matrix = zeros(Int64, num_categories, num_categories)
    for i in 1:num_categories
        for j in 1:i
            # Note: The following fails if i, j are not the same length.
            square_distance_matrix[i, j] = sum((codon_param_index_vec[i] .- codon_param_index_vec[j]) .^ 2)
            square_distance_matrix[j, i] = square_distance_matrix[i, j]  # Symmetric matrix
        end
    end
    return square_distance_matrix
end

function calculate_alloc_grid_and_theta(model::GeneralizedFUBARModel, ambient_samples::Vector{Vector{Float64}}, burnin::Int)
    """
    Samples the allocation grid from the ambient samples
    TODO: WHY DOES THIS TAKE SO MUCH TIME x_x
    """
    n_samples = size(ambient_samples, 1)
    n_sites = size(model.con_lik_matrix, 2)
    alloc_grid = zeros(Int64, size(model.con_lik_matrix))
    theta = zeros(Float64, model.n_categories)
    println()
    @showprogress for i = burnin+1:n_samples
        probability_vector = model.to_probability_vector(ambient_samples[i])
        theta .+= probability_vector # Accumulate the probability vector
        v = zeros(Float64, model.n_categories)
        for i = 1:n_sites
            for j = 1:model.n_categories
                v[j] = probability_vector[j] * model.con_lik_matrix[j, i]
            end
            samp = sample(1:model.n_categories, Weights(v))
            alloc_grid[samp, i] += 1 # Increment the allocation for the sampled category
        end
    end
    return alloc_grid, theta ./ (n_samples - burnin)
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
        # We should figure out a good prior dist for F(s)
        transition_function = s -> CodonMolecularEvolution.quintic_smooth_transition(s, -1, 1)
        # Define the masks for the suppression parameters
        hypothesis_masks = ones(Bool, (4, size(con_lik_matrix)[1]))
        hypothesis_masks[1, :] = [c[2] > 1 for c in codon_param_vec] # omega_1 > 1
        hypothesis_masks[2, :] = [c[3] > 1 for c in codon_param_vec] # omega_2 > 1
        hypothesis_masks[2, :] = [c[2] > c[3] for c in codon_param_vec] # omega_1 > omega_2
        hypothesis_masks[3, :] = [c[3] > c[2] for c in codon_param_vec] # omega_2 > omega_1

        square_distance_matrix = generate_square_l2_distance_matrix(codon_param_index_vec)
        kernel_stddev = 1.0 # example values idk what these should be xD
        suppression_stddev = 1.0
        model = SKBDIModel(
            [alphagrid, omegagrid, background_omega_grid],
            ["alpha", "omega_1", "omega_2", "background_omega"],
            hypothesis_masks,
            transition_function,
            log_con_lik_matrix,
            con_lik_matrix,
            codon_param_vec,
            codon_param_index_vec,
            ambient_sample -> hedwigs_ambient_to_parameter_transform(model,
                ambient_sample,
                kernel_stddev,
                suppression_stddev,
                square_distance_matrix,
                fast_cov_mat_hedwigs_kernel),
            1
        )
        # One further step of abstraction skbdi -> generalized fubar
        fubar_model = GeneralizedFUBARModel(model)
        ess_model = ESSModel(fubar_model.prior, fubar_model.log_likelihood)
        sample_time = @elapsed begin
            ambient_samples = AbstractMCMC.sample(ess_model, ESS(), iters, progress=true)
            alloc_grid, theta = calculate_alloc_grid_and_theta(fubar_model, ambient_samples, burnin)
        end
    end
    # Now we should have the same stuff that diffubar generates.
    df, plots_named_tuple = CodonMolecularEvolution.difFUBAR_tabulate_and_plot(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors, verbosity=verbosity, exports=exports)
    json = CodonMolecularEvolution.dNdS2JSON(CodonMolecularEvolution.difFUBAR2JSON(), (outpath=analysis_name, df=df, θ=theta, posterior_mat=alloc_grid ./ sum(alloc_grid, dims=1), categories=reduce(hcat, codon_param_vec)', tags=tags, tree=shallow_tree, LL=LL, timers=(total_time, fit_time, grid_time, sample_time), treestring=treestring, seqnames=seqnames, seqs=seqs, leaf_name_transform=leaf_name_transform, pos_thresh=pos_thresh, iters=iters, burnin=burnin, concentration=concentration, binarize=binarize, exports=exports2json))
    #Return df, (tuple of partial calculations needed to re-run tablulate), plots
    push!(plot_collection, plots_named_tuple)
    return df, (alloc_grid, ambient_samples, model, tag_colors), merge(plot_collection...)
end

analysis_name = "output/Ace2"
seqnames, seqs = read_fasta("test/data/Ace2_tiny/Ace2_tiny_tagged.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("test/data/Ace2_tiny/Ace2_tiny_tagged.nex")
df, results = skbdifFUBAR(seqnames, seqs, treestring, tags, analysis_name, tag_colors)
alloc_grid, ambient_samples, output_grid, _ = results
@save "output/alloc_grid.jld2" alloc_grid
@save "output/ambient_samples.jld2" ambient_samples
@save "output/grid.jld2" output_grid




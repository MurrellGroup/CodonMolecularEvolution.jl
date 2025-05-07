
#Takes in a vector of tags, and returns a function that will strip these from sequence names
function generate_tag_stripper(tags)
    function strip_tags_from_name(s::String)
        for t in tags
            s = replace(s, t => "")
        end
        return s
    end
    return strip_tags_from_name
end

#Rename this nonsense
"""
    import_colored_figtree_nexus_as_tagged_tree(fname; custom_labels=String[])

Takes a nexus file from FigTree, where branches have been colored. Replaces all color tags with group tags that can be used in the models. Can add custom labels too. Should consider an entire custom dictionary as well in future.

# Examples
```julia-repl
julia> treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/Ace2_no_background.nex")
("(((XM_027533928_Bos_indicus_x_Bos_taurus{G1}:0.097072,(XM_042974087_Panthera_tigris{G1}:0.038016,... more ...;", ["{G2}", "{G1}"], ["#ff0015", "#0011ff"])
```
!!! note
    `treestring` is truncated. [NEXUS tree file](https://raw.githubusercontent.com/MurrellGroup/CodonMolecularEvolution.jl/main/test/data/Ace2_no_background/Ace2_no_background.nex)
"""
function import_colored_figtree_nexus_as_tagged_tree(fname; custom_labels=String[])
    start_substr = "[&R] "
    lines = readlines(fname)
    treeline = lines[findfirst([occursin(start_substr, l) for l in lines])]
    st = findfirst(start_substr, treeline)[end]
    treestr = treeline[st+1:end]
    R = r"\[\&\!color\=\#\w*\]"
    color_tags = union([string(m.match) for m in eachmatch(R, treestr)])
    d = Dict{String,String}()

    if length(color_tags) == 0
        @warn "No color tags detected."
    else
        if length(custom_labels) == 0
            for i in 1:length(color_tags)
                d[color_tags[i]] = "{G$(i)}"
            end
        else
            for i in 1:length(color_tags)
                d[color_tags[i]] = custom_labels[i]
            end
        end
    end
    for k in keys(d)
        treestr = replace(treestr, k => d[k])
    end

    tag_colors = [k[10:end-1] for k in keys(d)]
    tags = collect(values(d))

    return treestr, tags, tag_colors
end

export import_colored_figtree_nexus_as_tagged_tree

########################################################
#Function for grouped Newick to DifFUBAR analysis format
########################################################

function import_hyphy_simulated_FASTA(file_name)
    # Import data generated using https://github.com/veg/hyphy-analyses/tree/master/SimulateMG94
    sequence_names = []
    sequences = []

    open(file_name, "r") do file
        current_sequence = ""

        for line in eachline(file)
            if length(line) == 0
                break
            else
                if line[1] == '>'
                    # Store the current sequence and reset for the next entry
                    if !isempty(current_sequence)
                        push!(sequences, current_sequence)
                    end
                    push!(sequence_names, line[2:end])
                    current_sequence = ""
                else
                    # Append the line to the current sequence
                    current_sequence *= line
                end
            end
        end

        # Store the last sequence after the loop ends (if any)
        if !isempty(current_sequence)
            push!(sequences, current_sequence)
        end
    end
    sequence_names = [string(seq_name) for seq_name in sequence_names]
    sequences = [string(seq) for seq in sequences]

    return sequence_names, sequences
end

export import_hyphy_simulated_FASTA

function replace_newick_tags(treestr)
    # This function replaces whatever is in {} in a newick tree to G1...Gn

    pattern = r"\{([^}]+)\}"
    unique_tags = Set{String}()
    for match in eachmatch(pattern, treestr)
        push!(unique_tags, match.match)
    end
    unique_tags = collect(unique_tags)
    num_unique_tags = length(unique_tags)
    group_tags = ["{G$index}" for (index, _) in enumerate(unique_tags)]

    if length(unique_tags) != length(group_tags)
        throw(ArgumentError("The number of unique tags must be equal to the number of group tags."))
    end

    tag_mapping = Dict{String,String}()

    for (old_tag, new_tag) in zip(unique_tags, group_tags)
        tag_mapping[old_tag] = "$new_tag"
    end

    for (old_tag, new_tag) in tag_mapping
        treestr = replace(treestr, old_tag => new_tag)
    end

    return treestr, group_tags, unique_tags
end


# Define a function to generate distinct hexadecimal color codes
function generate_hex_colors(num_colors)
    # Generates an arbitrary number of hex-color to match group_tags
    colors = []
    for i in 1:num_colors
        r = rand(0:255)
        g = rand(0:255)
        b = rand(0:255)
        hex_color = "#" * string(r, base=16, pad=2) * string(g, base=16, pad=2) * string(b, base=16, pad=2)
        push!(colors, hex_color)
    end
    colors_string = [string(color) for color in colors]
    return colors_string
end

"""
    import_grouped_label_tree(tree_file)

Takes a Newick tree file and return Newick tree, Newick tree with replaced tags, group tags, original tags, and randomly generated colours for each tag
"""
function import_grouped_label_tree(tree_file)
    tree = read_newick_tree(tree_file)
    treestring = newick(tree)
    treestring_group_labeled, group_tags, original_tags = replace_newick_tags(treestring)
    tag_colors = generate_hex_colors(length(original_tags))
    return treestring_group_labeled, treestring, group_tags, original_tags, tag_colors
end

export import_grouped_label_tree


function select_analysis_tags_from_newick_tree(tags, tag_colors, tag_pos)
    # If there are more than 2 tags in the newick tree, we split up the tags to the tags to be analyzed and the tags to be removed (placed as background)
    analysis_tags = tags[tag_pos]
    analysis_tag_colors = tag_colors[tag_pos]
    remove_tags = tags[filter(i -> i ∉ tag_pos, eachindex(tags))]
    return analysis_tags, analysis_tag_colors, remove_tags
end
export select_analysis_tags_from_newick_tree

function remove_tags_from_newick_tree(treestring, tags)
    # Cleans Newick tree by removing tags from tree before analysis step. Select subset of groups the rest is background. 
    new_tree = treestring
    for tag in tags
        new_tree = replace(new_tree, tag => "")
    end
    return new_tree
end
export remove_tags_from_newick_tree

export import_labeled_phylotree_newick
"""
    import_labeled_phylotree_newick(fname)

Import a tagged phylogeny from phylotree and return the treestring and tags.
"""
function import_labeled_phylotree_newick(fname)
    treestring = read(fname, String)
    R = r"(\{[^}]+\})"
    tags = union([string(m.match) for m in eachmatch(R, treestring)])
    return treestring, tags
end



#TODO: Make this first optimize a nuc model, then optimize a codon model starting from the nuc model estimates, fixing the GTR matrix
"""
    optimize_MG94_F3x4(seqnames, seqs, tree; leaf_name_transform=x -> x, genetic_code=MolecularEvolution.universal_code)

Optimizes the MG94+F3x4 model on a tree, given a set of sequences and a tree. Returns the optimized tree, alpha, beta, nuc_matrix, F3x4, and eq_freqs.
The leaf_name_transform kwarg can be used to transform the leaf names in the tree to match the seqnames.
"""
function optimize_MG94_F3x4(seqnames, seqs, tree; leaf_name_transform=x -> x, genetic_code=MolecularEvolution.universal_code)

    #Count F3x4 frequencies from the seqs, and estimate codon freqs from this
    f3x4 = MolecularEvolution.count_F3x4(seqs)
    eq_freqs = MolecularEvolution.F3x4_eq_freqs(f3x4)

    #Set up a codon partition (will default to Universal genetic code)
    eq_partition = CodonPartition(Int64(length(seqs[1]) / 3), code=genetic_code)
    eq_partition.state .= eq_freqs
    initial_partition = LazyPartition{CodonPartition}()
    populate_tree!(tree, initial_partition, seqnames, seqs, leaf_name_transform=leaf_name_transform)
    lazyprep!(tree, [eq_partition])

    #We'll use the empirical F3x4 freqs, fixed MG94 alpha=1, and optimize the nuc parameters and MG94 beta
    #Note: the nuc rates are confounded with alpha
    initial_params = (
        rates=positive(ones(6)), #rates must be non-negative
        beta=positive(1.0)
    )
    flat_initial_params, unflatten = value_flatten(initial_params) #See ParameterHandling.jl docs
    num_params = length(flat_initial_params)

    function build_model_vec(p; F3x4=f3x4, alpha=1.0)
        #Need to pass through genetic code here!
        #If you run into numerical issues with DiagonalizedCTMC, switch to GeneralCTMC instead
        return GeneralCTMC(MolecularEvolution.MG94_F3x4(alpha, p.beta, reversibleQ(p.rates, ones(4)), F3x4, genetic_code=genetic_code))
    end

    function objective(params::NamedTuple; tree=tree, eq_freqs=eq_freqs)
        return -log_likelihood!(tree, build_model_vec(params))
    end

    opt = Opt(:LN_BOBYQA, num_params)
    min_objective!(opt, (x, y) -> (objective ∘ unflatten)(x))
    lower_bounds!(opt, [-5.0 for i in 1:num_params])
    upper_bounds!(opt, [5.0 for i in 1:num_params])
    xtol_rel!(opt, 1e-12)
    _, mini, _ = NLopt.optimize(opt, flat_initial_params)

    final_params = unflatten(mini)

    #tree, alpha, beta, nuc_matrix, F3x4, eq_freqs
    return tree, 1.0, final_params.beta, reversibleQ(final_params.rates, ones(4)), f3x4, eq_freqs
end

export optimize_MG94_F3x4

function optimize_nuc_mu(seqnames, seqs, tree; leaf_name_transform=x -> x, genetic_code=MolecularEvolution.universal_code, optimize_branch_lengths = false)
    #Optimize mu params using a nuc model, mu denotes the nucleotide mutational biases
    nuc_pi = char_proportions(seqs, MolecularEvolution.nucstring)

    eq_partition = NucleotidePartition(length(seqs[1]))
    eq_partition.state .= nuc_pi
    if optimize_branch_lengths == true
        populate_tree!(tree, eq_partition, seqnames, seqs, leaf_name_transform=leaf_name_transform)
    else
        initial_partition = LazyPartition{NucleotidePartition}()
        populate_tree!(tree, initial_partition, seqnames, seqs, leaf_name_transform=leaf_name_transform)
        lazyprep!(tree, [eq_partition])
    end

    initial_params = positive(ones(6)) #rates must be non-negative

    flat_initial_params, unflatten = value_flatten(initial_params) #See ParameterHandling.jl docs
    num_params = length(flat_initial_params)

    function build_model_vec(p; nuc_pi=nuc_pi)
        return GeneralCTMC(reversibleQ(p, nuc_pi))
    end

    function objective(params; tree=tree)
        return -log_likelihood!(tree, build_model_vec(params))
    end

    opt = Opt(:LN_BOBYQA, num_params)
    min_objective!(opt, (x, y) -> (objective ∘ unflatten)(x))
    lower_bounds!(opt, [-5.0 for i in 1:num_params])
    upper_bounds!(opt, [5.0 for i in 1:num_params])
    xtol_rel!(opt, 1e-12)
    _, mini, _ = NLopt.optimize(opt, flat_initial_params)

    final_params = unflatten(mini)
    #tree, nuc mu rates, nuc pi rates
    return tree, final_params, nuc_pi
end

function optimize_codon_alpha_and_beta(seqnames, seqs, tree, GTRmat; leaf_name_transform=x -> x, genetic_code=MolecularEvolution.universal_code, verbosity = 1)
    #Now we optimize alpha and beta rates using a codon model
    #Count F3x4 frequencies from the seqs, and estimate codon freqs from this
    f3x4 = MolecularEvolution.count_F3x4(seqs)
    eq_freqs = MolecularEvolution.F3x4_eq_freqs(f3x4)

    #Set up a codon partition (will default to Universal genetic code)
    eq_partition = CodonPartition(Int64(length(seqs[1]) / 3), code=genetic_code)
    eq_partition.state .= eq_freqs
    initial_partition = LazyPartition{CodonPartition}()
    populate_tree!(tree, initial_partition, seqnames, seqs, leaf_name_transform=leaf_name_transform)
    lazyprep!(tree, [eq_partition])

    function build_model_vec(alpha, beta; F3x4=f3x4)
        #Need to pass through genetic code here!
        #If you run into numerical issues with DiagonalizedCTMC, switch to GeneralCTMC instead
        return robust_CTMC(MolecularEvolution.MG94_F3x4(alpha, beta, GTRmat, F3x4, genetic_code=genetic_code))
    end

    function objective(alpha=1, beta=1; tree=tree, eq_freqs=eq_freqs)
        return log_likelihood!(tree, build_model_vec(alpha, beta))
    end
    
    #High tolerance estimates of alpha then beta
    alpha, beta = 1.0, 1.0
    num_1d_optims = 3 #We'll do this many 1D optimizations, first for alpha, then beta, then alpha again and so forth
    initial_param = 1.0 #rate must be non-negative
    lower_bound = 0.0
    upper_bound = 5.0
    high_tol = 1e-4
    for i in 1:num_1d_optims-1
        isodd = i % 2 == 1
        if isodd
            alpha = brents_method_minimize(x -> -objective(x, beta), lower_bound, upper_bound, identity, high_tol)
            #@show Int(ceil(log(high_tol / (upper_bound - lower_bound)) / log(MolecularEvolution.invphi)))
            #alpha = golden_section_maximize(x -> objective(x, beta), lower_bound, upper_bound, identity, high_tol)
            #@show alpha1, alpha
        else
            beta = brents_method_minimize(x -> -objective(alpha, x), lower_bound, upper_bound, identity, high_tol)
            #beta = golden_section_maximize(x -> objective(alpha, x), lower_bound, upper_bound, identity, high_tol)
            #@show beta1, beta
        end
    end

    #Low tolerance final estimate of alpha
    """
    the parabola thing is: 
    take 3 values of the parameter, 
    and eval the LL. fit the parabola, 
    and take the max of the parabola, 
    then swap out the worst LL in the previous set for the new max, 
    and repeat with those 3 values 
    // Benjamin Murrell in a teams convo
    """
    low_tol = 1e-7
    ε = 1e-2
    prev_alpha = Inf
    maxiters = 20 #We should be quite close to the final alpha
    iters = 0
    xvec = [alpha - ε, alpha, alpha + ε]
    yvec = map(x -> objective(x, beta), xvec)
    max_obj = 0.0
    while abs(alpha - prev_alpha) > low_tol && iters < maxiters
        prev_alpha = alpha
        if iters > 0 #We use this order to avoid one potential objective call
            wipe_index = argmin(yvec)
            max_obj = objective(alpha, beta)
            xvec[wipe_index] = alpha
            yvec[wipe_index] = max_obj
            perm = sortperm(xvec)
            xvec = xvec[perm]
            yvec = yvec[perm]  
        end
        alpha = MolecularEvolution.new_max(xvec, yvec)
        iters += 1
    end
    
    verbosity > 0 && println("Optimized single α,β LL=$(max_obj) with α=$(alpha) and β=$(beta).")
    #tree, alpha, beta, F3x4, eq_freqs
    return tree, alpha, beta, f3x4, eq_freqs
end

function rescale_branchlengths!(tree, scale)
    for node in getnodelist(tree)
        node.branchlength *= scale
    end
end

##########################################
#Possibly migrate to MolecularEvolution.jl
##########################################

#This should be in the main repo
function LDA_gibbs_track_allocation_vec(conditionals::Array{Float64,2}, alpha::Float64; iters=10000)
    grid_size, num_sites = size(conditionals)
    #This should instead track an integer allocation vec per MCMC iteration - will be smaller
    alloc_grid = zeros(Int64, (grid_size, num_sites))
    φ = zeros(grid_size)
    θ = ones(grid_size) ./ grid_size
    v = zeros(grid_size)
    θsum = zeros(grid_size)
    burnin::Int = div(iters, 5)
    for iter = 1:iters
        φ .= alpha
        for i = 1:num_sites
            @simd for j = 1:grid_size
                @inbounds v[j] = θ[j] * conditionals[j, i]
            end
            samp = sample(1:grid_size, Weights(v))
            φ[samp] += 1.0
            if iter > burnin
                alloc_grid[samp, i] += 1
            end
        end
        θ .= rand(Dirichlet(φ))
        if iter > burnin
            #NOTE: test without simd and inbounds macros
            @simd for j = 1:grid_size
                @inbounds θsum[j] += θ[j]
            end
        end
    end
    return alloc_grid, θsum ./ sum(θsum)
end

#Tries to init a DiagonalizedCTMC and falls back to GeneralCTMC if it fails
function robust_CTMC(Q; kwargs...)
    try
        return DiagonalizedCTMC(Q; kwargs...)
    catch
        print(".")
        return GeneralCTMC(Q; kwargs...)
    end
end

#Sets up model memoization - could probably do this a bit more generally
function MG94_cacher(code)
    d = Dict{Any,DiagonalizedCTMC}()
    function cached_model(args...; genetic_code=code)
        if !haskey(d, args)
            d[args] = DiagonalizedCTMC(MolecularEvolution.MG94_F3x4(args..., genetic_code=genetic_code))
        end
        return d[args]
    end
    return cached_model
end

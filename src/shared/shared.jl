
#Takes in a vector of tags, and returns a function that will strip these from sequence names
function generate_tag_stripper(tags)
    function strip_tags_from_name(s::String)
        for t in tags
            s = replace(s, t=>"")
        end
        return s
    end
    return strip_tags_from_name
end

#Rename this nonsense
"""
    scrape_figtree_colored_nexus(fname; custom_labels = String[])

Takes a nexus file from FigTree, where branches have been colored. Replaces all color tags with group tags that can be used in the models. Can add custom labels too. Should consider an entire custom dictionary as well in future.
"""
function import_colored_figtree_nexus_as_tagged_tree(fname; custom_labels = String[])
    start_substr = "[&R] ";
    lines = readlines(fname);
    treeline = lines[findfirst([occursin(start_substr,l) for l in lines])];
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
        treestr = replace(treestr,k=>d[k])
    end

    tag_colors = [k[10:end-1] for k in keys(d)]
    tags = collect(values(d))

    return treestr, tags, tag_colors
end

export import_colored_figtree_nexus_as_tagged_tree


#TODO: Make this first optimize a nuc model, then optimize a codon model starting from the nuc model estimates, fixing the GTR matrix
"""
    optimize_MG94_F3x4(seqnames, seqs, treestring; leaf_name_transform = x -> x)

Optimizes the MG94+F3x4 model on a tree, given a set of sequences and a tree. Returns the optimized tree, alpha, beta, nuc_matrix, F3x4, and eq_freqs.
The leaf_name_transform kwarg can be used to transform the leaf names in the tree to match the seqnames.
"""
function optimize_MG94_F3x4(seqnames, seqs, tree; leaf_name_transform = x -> x, genetic_code = MolecularEvolution.universal_code)
    
    #Count F3x4 frequencies from the seqs, and estimate codon freqs from this
    f3x4 = MolecularEvolution.count_F3x4(seqs)
    eq_freqs = MolecularEvolution.F3x4_eq_freqs(f3x4)

    #Set up a codon partition (will default to Universal genetic code)
    initial_partition = CodonPartition(Int64(length(seqs[1])/3), code = genetic_code)
    initial_partition.state .= eq_freqs
    populate_tree!(tree,initial_partition,seqnames,seqs, leaf_name_transform = leaf_name_transform)

    #We'll use the empirical F3x4 freqs, fixed MG94 alpha=1, and optimize the nuc parameters and MG94 beta
    #Note: the nuc rates are confounded with alpha
    initial_params = (
            rates = positive(ones(6)), #rates must be non-negative
            beta = positive(1.0)
    )
    flat_initial_params, unflatten = value_flatten(initial_params) #See ParameterHandling.jl docs
    num_params = length(flat_initial_params)

    function build_model_vec(p; F3x4 = f3x4, alpha = 1.0)
        #Need to pass through genetic code here!
        #If you run into numerical issues with DiagonalizedCTMC, switch to GeneralCTMC instead
        return GeneralCTMC(MolecularEvolution.MG94_F3x4(alpha, p.beta, reversibleQ(p.rates,ones(4)), F3x4, genetic_code = genetic_code))
    end

    function objective(params::NamedTuple; tree = tree, eq_freqs = eq_freqs)
        return -log_likelihood!(tree,build_model_vec(params))
    end

    opt = Opt(:LN_BOBYQA, num_params)
    min_objective!(opt, (x,y) -> (objective ∘ unflatten)(x))
    lower_bounds!(opt, [-5.0 for i in 1:num_params])
    upper_bounds!(opt, [5.0 for i in 1:num_params])
    xtol_rel!(opt, 1e-12)
    _,mini,_ = NLopt.optimize(opt, flat_initial_params)

    final_params = unflatten(mini)
    
    #tree, alpha, beta, nuc_matrix, F3x4, eq_freqs
    return tree,1.0,final_params.beta,reversibleQ(final_params.rates,ones(4)),f3x4,eq_freqs 
end




##########################################
#Possibly migrate to MolecularEvolution.jl
##########################################

#This should be in the main repo
function LDA_gibbs_track_allocation_vec(conditionals::Array{Float64,2}, alpha::Float64; iters = 10000)
    grid_size, num_sites = size(conditionals)
    #This should instead track an integer allocation vec per MCMC iteration - will be smaller
    alloc_grid = zeros(Int64,(grid_size,num_sites))
    φ = zeros(grid_size)
    θ = ones(grid_size) ./ grid_size
    v = zeros(grid_size)
    θsum = zeros(grid_size)
    burnin::Int = div(iters, 5)
    for iter=1:iters
        φ .= alpha
        for i=1:num_sites
                @simd for j=1:grid_size
                @inbounds v[j] = θ[j]*conditionals[j,i]
                end
            samp = sample(1:grid_size,Weights(v))
            φ[samp] += 1.0
            if iter > burnin
                alloc_grid[samp,i] += 1
            end
        end
        θ .= rand(Dirichlet(φ))
        if iter > burnin
            #NOTE: test without simd and inbounds macros
            @simd for j=1:grid_size
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
    function cached_model(args...; genetic_code = code)
        if !haskey(d,args)
            d[args] = DiagonalizedCTMC(MolecularEvolution.MG94_F3x4(args..., genetic_code = genetic_code))
        end
        return d[args]
    end
    return cached_model
end

#Modified to allow an option renaming function that can eg. strip tags from the tree when matching seq data
function MolecularEvolution.populate_tree!(
    tree::FelNode,
    starting_message::Vector{<:Partition},
    names,
    data;
    init_all_messages = true,
    tolerate_missing = 1, #0 = error if missing; 1 = warn and set to missing data; 2 = set to missing data
    leaf_name_transform = x -> x
)
    if init_all_messages
        internal_message_init!(tree, starting_message)
    else
        tree.parent_message = deepcopy(starting_message)
    end
    name_dic = Dict(zip(names, 1:length(names)))
    for n in getleaflist(tree)
        if haskey(name_dic, leaf_name_transform(n.name))
            MolecularEvolution.populate_message!(n.message, data[name_dic[leaf_name_transform(n.name)]])
        else
            warn_str = n.name * " on tree but not found in names."
            if tolerate_missing == 0
                @error warn_str
            end
            if tolerate_missing == 1
                @warn warn_str
            end
            uninformative_message!(n.message)
        end
    end
end

function MolecularEvolution.populate_tree!(
    tree::FelNode,
    starting_partition::Partition,
    names,
    data;
    init_all_messages = true,
    tolerate_missing = 1,
    leaf_name_transform = x -> x
)
    MolecularEvolution.populate_tree!(
        tree,
        [starting_partition],
        names,
        data,
        init_all_messages = init_all_messages,
        tolerate_missing = tolerate_missing,
        leaf_name_transform = leaf_name_transform
    )
end


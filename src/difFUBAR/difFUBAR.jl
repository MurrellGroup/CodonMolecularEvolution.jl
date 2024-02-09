#Things to think about:
#1 - A standardized interface for all of our user-facing models.

#Thids to do:
#1 - If the only node without a tag is the root node, then this should still count as the "no-background" case
#2 - Add code that automatically infers internal node tags based on:
#  A - simple discrete models
#  B - simple continuous models
#3 - Add more options for initial model fit
#4 - Get subtree-caching working
#5 - Add tree import where the user specifies FG1, and everything else is FG2 (ie. no background)


#Here, untagged comes after the tags
function model_ind(str::String, tags::Vector{String})
    ind = length(tags)+1
    for (i,t) in enumerate(tags)
        if occursin(t,str)
            ind = i
        end
    end
    return ind
end

function add_to_each_element(vec_of_vec, elems)
    return [vcat(v,[e]) for v in vec_of_vec for e in elems]
end

#This is the function that assigns models to branches
#Sometimes there will be the same number of tags as omegas, but sometimes there will be one more omega.
#Depending on whether there are any background branches (UNTESTED WITH NO BACKGROUND BRANCHES)
function N_Omegas_model_func(cached_model, tags,omega_vec,alpha,nuc_mat,F3x4, code)
    models = [cached_model(alpha,alpha*o,nuc_mat,F3x4, genetic_code = code) for o in omega_vec];
    return n::FelNode -> [models[model_ind(n.name,tags)]]
end

function Omega_model_func(cached_model,omega,alpha,nuc_mat,F3x4, code)
    model = cached_model(alpha,alpha*omega,nuc_mat,F3x4, genetic_code = code);
    return n::FelNode -> [model]
end

function collapse_counts(param_vec,count_vec; cases = nothing)
    if isnothing(cases)
        cases = sort(union(param_vec))
    end
    d = Dict(zip(cases,1:length(cases)))
    storage = zeros(Int64,length(cases))
    for i in 1:length(count_vec)
        storage[d[param_vec[i]]] += count_vec[i]
    end
    return storage ./ sum(storage)
end

function FUBAR_violin_plot(sites,group1_volumes,omegagrid;
    color = "black", tag = "", alpha = 0.6, 
    x_label = "Parameter",y_label = "Codon Sites",
    v_offset = 0.0, legend_ncol = 3,
    vertical_ind = findfirst(omegagrid .>= 1.0), 
    plot_legend = true)
    ypos = [-i*0.5 for i in 1:length(sites)]
    if !isnothing(vertical_ind)
        bar([vertical_ind],[2 + maximum(ypos)-minimum(ypos)], bottom = [minimum(ypos) - 1] , color = "grey", alpha = 0.05)
    end
    for i in 1:length(sites)
        center_line = ypos[i]
        bar((1:length(omegagrid)),group1_volumes[i], bottom = v_offset .+ center_line .+ (-0.5 .* group1_volumes[i]), color = color, alpha = alpha)    
    end    
    bar([1],[1], bottom = [1000] , color = color, alpha = alpha, label = tag)
    yticks(ypos,sites)
    xticks((1:length(omegagrid)),omegagrid, rotation = 90);
    ylabel(y_label)
    xlabel(x_label)

    ylim(minimum(ypos)-0.5, maximum(ypos)+0.5)
    if plot_legend
        legend(loc="upper center", bbox_to_anchor=(0.5, 1.01),
              ncol=legend_ncol, fancybox=true, shadow=true)
    end
    tight_layout()
end


"""
    getpuresubclades(node::FelNode, tags::Vector{String}, pure_subclades=FelNode[])

- Should usually be called on the root of the tree. Traverses the tree recursively with a depth-first search to find roots of pure subclades, presuming that nodenames have been trailed with tags.
- To just get the pure subclades, one can run `pure_subclades, _, _ = getpuresubclades(tree, tags)`.

# Arguments
- `node`: The root of the search.
- `tags`: A vector of tags.
- `pure_subclades`: A vector of pure subclades. Defaults to an empty vector.

# Returns
- A tuple containing the vector of pure subclades, a boolean indicating whether the node is pure, and the index of the node's tag.
"""
function getpuresubclades(node::FelNode, tags::Vector{String}, pure_subclades=FelNode[])
    # Get the index of the node's tag
    tag_ind_of_node = model_ind(node.name, tags)

    # If the node is a leaf, it's pure
    if isleafnode(node)
        return pure_subclades, true, tag_ind_of_node
    end

    children_are_pure = Vector{Bool}()
    children_tag_inds = Vector{Int64}()

    for child in node.children
        pure_subclades, child_is_pure, tag_ind = getpuresubclades(child, tags, pure_subclades)
        push!(children_are_pure, child_is_pure)
        push!(children_tag_inds, tag_ind)
    end

    # Get the index of the node's first child's tag
    tag_ind_of_first_child = first(children_tag_inds)

    # This is the case where the subclade starting at node is pure
    if all(children_are_pure) && all(x == tag_ind_of_first_child for x in children_tag_inds)
        if tag_ind_of_node != tag_ind_of_first_child
            # The purity is broken at this node
            push!(pure_subclades, node)
            return pure_subclades, false, tag_ind_of_node
        end
        # The purity is not broken at this node
        return pure_subclades, true, tag_ind_of_node
    end

    # This is the case where some child has mixed tags or the children are pure with regards to different tags
    for (child_is_pure, child) in zip(children_are_pure, node.children)
        if !child_is_pure || isleafnode(child)
            # We don't want to push leaves into pure_subclades
            continue
        end
        push!(pure_subclades, child)
    end
    return pure_subclades, false, tag_ind_of_node
end

"""


"""
function difFUBAR_init(outpath_and_file_prefix, treestring, tags, tag_colors; verbosity = 1, exports = true, strip_tags_from_name = generate_tag_stripper(tags))
    
    #Create the export directory, if required
    analysis_name = outpath_and_file_prefix
    splt = splitpath(analysis_name)[1:end-1]
    if length(splt) > 0
        exports && mkpath(joinpath(splt))
    end

    tree = gettreefromnewick(treestring, FelNode)

    #Need to think of consequences of eg. binarizing when there are tags.
    #We'll need the non-zero branch lengths to inherit the names/tags, for example.
    #REQUIRED TEST: CHECK NODE NAMES AFTER BINARIZATION
    MolecularEvolution.binarize!(tree) #Check if this is required for trees with ternary branching?
    #Make this optional, but done by default
    MolecularEvolution.ladderize!(tree)

    #data validation checks:
    #seq lengths must all be multiples of there
    #seqs must not contain stop codons, unless "handle_stop_codons" or "warn_stop_codons" is set
    #seqs must not contain ambiguous characters, unless "handle_ambiguous" or "warn_ambiguous" is set
    #seqs must not contain non-ATGC characters, unless "handle_nonATGC" or "warn_nonATGC" is set
    #tree names must match sequence names (after removing tags)

    #Export tagged input tree to file?
    verbosity > 0 && println("Step 1: Initialization. If exports = true, tree showing the assignment of branches to groups/colors will be exported to: "*analysis_name*"_tagged_input_tree.svg.")


    p = sortperm(tags)
    tags,tag_colors = tags[p],tag_colors[p]
    push!(tag_colors,"black") #ASSUMPTION: background color is always black

    if exports
        #Replace with Phylo.jl based plot?
        color_dict = Dict(zip(getnodelist(tree),[tag_colors[model_ind(n.name,tags)] for n in getnodelist(tree)]));
        label_dict = Dict(zip(getnodelist(tree),[strip_tags_from_name(n.name) for n in getnodelist(tree)]))
        img = tree_draw(tree, canvas_height = (3+length(getleaflist(tree))/5)cm,
            draw_labels = true, dot_color_dict = color_dict,
            line_color_dict = color_dict, line_width = 0.3, min_dot_size = 0.01,
        nodelabel_dict = label_dict)
        img |> SVG(analysis_name*"_tagged_input_tree.svg",15cm, (3+length(getleaflist(tree))/5)cm)
    end

    #Tags and tag colors are now ordered, and tag_colors includes the untagged category
    return tree, tags, tag_colors, analysis_name
end

function difFUBAR_global_fit(seqnames, seqs, tree, leaf_name_transform, code; verbosity = 1)

    verbosity > 0 && println("Step 2: Optimizing global codon model parameters.")

    tree,alpha,beta,GTRmat,F3x4_freqs,eq_freqs = optimize_MG94_F3x4(seqnames, seqs, tree, leaf_name_transform = leaf_name_transform, genetic_code = code)

    ######
    #optionally polish branch lengths and topology
    ######

    return tree, alpha,beta,GTRmat,F3x4_freqs,eq_freqs
end

#foreground_grid and background_grid control the number of categories below 1.0
function difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code; verbosity = 1, foreground_grid = 6, background_grid = 4)

    cached_model = MG94_cacher(code)

    #This is the function that assigns models to branches
    #Sometimes there will be the same number of tags as omegas, but sometimes there will be one more omega.
    #Depending on whether there are any background branches (UNTESTED WITH NO BACKGROUND BRANCHES)
    function N_Omegas_model_func(tags,omega_vec,alpha,nuc_mat,F3x4, code)
        models = [cached_model(alpha,alpha*o,nuc_mat,F3x4, genetic_code = code) for o in omega_vec];
        return n::FelNode -> [models[model_ind(n.name,tags)]]
    end

    #Defines the grid used for inference.
    function gridsetup(lb, ub, num_below_one, trin, tr)
        step = (trin(1.0) - trin(lb))/num_below_one
        return tr.(trin(lb):step:trin(ub))
    end
    tr(x) = 10^x-0.05
    trinv(x) =  log10(x+0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr); 
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree) if !isroot(n)]) > num_groups
    tensor_dims = 1+num_groups+is_background;

    function add_to_each_element(vec_of_vec, elems)
        return [vcat(v,[e]) for v in vec_of_vec for e in elems]
    end
    
    codon_param_vec = [[a] for a in alphagrid]
    param_kinds = ["Alpha"]
    for g in 1:num_groups
        push!(param_kinds, "OmegaG$(g)")
        codon_param_vec = add_to_each_element(codon_param_vec,omegagrid)
    end
    if is_background
        push!(param_kinds, "OmegaBackground")
        codon_param_vec = add_to_each_element(codon_param_vec,background_omega_grid)
    end
    codon_param_vec;
    
    num_sites = tree.message[1].sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l,num_sites);

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")

    for (row_ind,cp) in enumerate(codon_param_vec)
        alpha = cp[1]
        omegas = cp[2:end]
        tagged_models = N_Omegas_model_func(tags,omegas,alpha,GTRmat,F3x4_freqs, code)

        felsenstein!(tree,tagged_models)
        #This combine!() is needed because the current site_LLs function applies to a partition
        #And after a felsenstein pass, you don't have the eq freqs factored in.
        #We could make a version of log_likelihood() that returns the partitions instead of just the sum
        combine!.(tree.message,tree.parent_message)

        log_con_lik_matrix[row_ind,:] .= MolecularEvolution.site_LLs(tree.message[1]) #Check that these grab the scaling constants as well!
        verbosity > 0 && if mod(row_ind,500)==1
            print(round(100*row_ind/length(codon_param_vec)),"% ")
            flush(stdout)
        end
    end
    verbosity > 0 && println()

    con_lik_matrix = zeros(size(log_con_lik_matrix));
    site_scalers = maximum(log_con_lik_matrix, dims = 1);
    for i in 1:num_sites
        con_lik_matrix[:,i] .= exp.(log_con_lik_matrix[:,i] .- site_scalers[i])
    end

    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end

function difFUBAR_grid_parallel(tree, tags, GTRmat, F3x4_freqs, code; verbosity = 1, foreground_grid = 6, background_grid = 4)

    cached_models = [MG94_cacher(code) for _ = 1:Threads.nthreads()]
    trees = [tree, [deepcopy(tree) for _ = 1:(Threads.nthreads() - 1)]...]

    #Defines the grid used for inference.
    function gridsetup(lb, ub, num_below_one, trin, tr)
        step = (trin(1.0) - trin(lb))/num_below_one
        return tr.(trin(lb):step:trin(ub))
    end
    tr(x) = 10^x-0.05
    trinv(x) =  log10(x+0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr); 
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree) if !isroot(n)]) > num_groups
    tensor_dims = 1+num_groups+is_background;
    
    codon_param_vec = [[a] for a in alphagrid]
    param_kinds = ["Alpha"]
    for g in 1:num_groups
        push!(param_kinds, "OmegaG$(g)")
        codon_param_vec = add_to_each_element(codon_param_vec,omegagrid)
    end
    if is_background
        push!(param_kinds, "OmegaBackground")
        codon_param_vec = add_to_each_element(codon_param_vec,background_omega_grid)
    end
    codon_param_vec;
    
    num_sites = tree.message[1].sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l,num_sites);

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")

    Threads.@threads for (row_ind,cp) in collect(enumerate(codon_param_vec))
        id = Threads.threadid()
        tree_local = trees[id]
        cached_model = cached_models[id]
        alpha = cp[1]
        omegas = cp[2:end]
        tagged_models = N_Omegas_model_func(cached_model, tags,omegas,alpha,GTRmat,F3x4_freqs, code)

        felsenstein!(tree_local,tagged_models)
        #This combine!() is needed because the current site_LLs function applies to a partition
        #And after a felsenstein pass, you don't have the eq freqs factored in.
        #We could make a version of log_likelihood() that returns the partitions instead of just the sum
        combine!.(tree_local.message,tree_local.parent_message)
        # kanske inte ens behöver lock i samband med pre-allocation
        log_con_lik_matrix[row_ind,:] .= MolecularEvolution.site_LLs(tree_local.message[1]) #Check that these grab the scaling constants as well!
        #verbosity > 0 && if mod(row_ind,500)==1
        #    print(round(100*row_ind/length(codon_param_vec)),"% ")
        #    flush(stdout)
        #end
    end
    verbosity > 0 && println()

    con_lik_matrix = zeros(size(log_con_lik_matrix));
    site_scalers = maximum(log_con_lik_matrix, dims = 1);
    for i in 1:num_sites
        con_lik_matrix[:,i] .= exp.(log_con_lik_matrix[:,i] .- site_scalers[i])
    end

    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end

function difFUBAR_grid_chunks(tree, tags, GTRmat, F3x4_freqs, code; verbosity = 1, foreground_grid = 6, background_grid = 4)

    cached_models = [MG94_cacher(code) for _ = 1:Threads.nthreads()]
    trees = [tree, [deepcopy(tree) for _ = 1:(Threads.nthreads() - 1)]...]

    #Defines the grid used for inference.
    function gridsetup(lb, ub, num_below_one, trin, tr)
        step = (trin(1.0) - trin(lb))/num_below_one
        return tr.(trin(lb):step:trin(ub))
    end
    tr(x) = 10^x-0.05
    trinv(x) =  log10(x+0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr); 
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree) if !isroot(n)]) > num_groups
    tensor_dims = 1+num_groups+is_background;
    
    codon_param_vec = [[a] for a in alphagrid]
    param_kinds = ["Alpha"]
    for g in 1:num_groups
        push!(param_kinds, "OmegaG$(g)")
        codon_param_vec = add_to_each_element(codon_param_vec,omegagrid)
    end
    if is_background
        push!(param_kinds, "OmegaBackground")
        codon_param_vec = add_to_each_element(codon_param_vec,background_omega_grid)
    end
    codon_param_vec;
    
    num_sites = tree.message[1].sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l,num_sites);

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")

    function do_subgrid(tree, cached_model, cpv_chunk)
        # Note that cpv_chunk is already enumerated
        for (row_ind,cp) in cpv_chunk
            alpha = cp[1]
            omegas = cp[2:end]
            tagged_models = N_Omegas_model_func(cached_model, tags,omegas,alpha,GTRmat,F3x4_freqs, code)

            felsenstein!(tree,tagged_models)
            #This combine!() is needed because the current site_LLs function applies to a partition
            #And after a felsenstein pass, you don't have the eq freqs factored in.
            #We could make a version of log_likelihood() that returns the partitions instead of just the sum
            combine!.(tree.message,tree.parent_message)
            log_con_lik_matrix[row_ind,:] .= MolecularEvolution.site_LLs(tree.message[1]) #Check that these grab the scaling constants as well!
            #verbosity > 0 && if mod(row_ind,500)==1
            #    print(round(100*row_ind/length(codon_param_vec)),"% ")
            #    flush(stdout)
            #end
        end
    end

    cpv_chunks = Iterators.partition(enumerate(codon_param_vec), max(1, ceil(Int, length(codon_param_vec) / Threads.nthreads())))
    tasks = []
    for (i, cpv_chunk) in enumerate(cpv_chunks)
        # Spawn the task and add it to the array
        task = Threads.@spawn do_subgrid(trees[i], cached_models[i], cpv_chunk)
        push!(tasks, task)
    end

    # Wait for all tasks to finish
    foreach(wait, tasks)

    #@show length(first(cpv_chunks))
    # cp = codon_param_vec[1]
    # @show trees[1].message[1].state[1]
    # @show trees[2].message[1].state[1]
    # keys1 = map(x -> (x[1]), collect(keys(cached_models[1](nothing, access_dict=true))))
    # @show unique(keys1)
    # println()
    # keys2 = map(x -> (x[1]), collect(keys(cached_models[2](nothing, access_dict=true))))
    # @show unique(keys2)
    # @show intersect(keys1, keys2)
    # @show codon_param_vec[max(1, length(codon_param_vec) ÷ Threads.nthreads()) + 1][1]

    verbosity > 0 && println()

    con_lik_matrix = zeros(size(log_con_lik_matrix));
    site_scalers = maximum(log_con_lik_matrix, dims = 1);
    for i in 1:num_sites
        con_lik_matrix[:,i] .= exp.(log_con_lik_matrix[:,i] .- site_scalers[i])
    end

    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end

function difFUBAR_grid_tree_surgery(tree, tags, GTRmat, F3x4_freqs, code; verbosity = 1, foreground_grid = 6, background_grid = 4)

    MolecularEvolution.set_node_indices!(tree)
    cached_model = MG94_cacher(code)

    #Defines the grid used for inference.
    function gridsetup(lb, ub, num_below_one, trin, tr)
        step = (trin(1.0) - trin(lb))/num_below_one
        return tr.(trin(lb):step:trin(ub))
    end
    tr(x) = 10^x-0.05
    trinv(x) =  log10(x+0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr); 
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree) if !isroot(n)]) > num_groups
    @show is_background
    tensor_dims = 1+num_groups+is_background;
    
    codon_param_vec = [[a] for a in alphagrid]
    param_kinds = ["Alpha"]
    for g in 1:num_groups
        push!(param_kinds, "OmegaG$(g)")
        codon_param_vec = add_to_each_element(codon_param_vec,omegagrid)
    end
    if is_background
        push!(param_kinds, "OmegaBackground")
        codon_param_vec = add_to_each_element(codon_param_vec,background_omega_grid)
    end
    codon_param_vec;

    num_sites = tree.message[1].sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l,num_sites);

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")

    pure_subclades, _, _ = getpuresubclades(tree, tags)

    if length(pure_subclades) > 0
        alpha_and_single_omega_grids = Dict()
        alphagrid_vectorized = [[a] for a in alphagrid]
        alpha_and_single_omega_grids["Omega"] = add_to_each_element(alphagrid_vectorized,omegagrid)
        if is_background
            alpha_and_single_omega_grids["OmegaBackground"] = add_to_each_element(alphagrid_vectorized,background_omega_grid)
        end
    end

    cached_messages = Dict()
    cached_tag_inds = Dict()
    model_time = 0
    copy_time = 0
    
    for x in pure_subclades
        @time begin
            tag_ind_below = model_ind(x.children[1].name, tags)
            nodeindex = x.nodeindex
            cached_tag_inds[nodeindex] = tag_ind_below
            if tag_ind_below <= num_groups
                alpha_and_single_omega_grid = alpha_and_single_omega_grids["Omega"]
            else
                alpha_and_single_omega_grid = alpha_and_single_omega_grids["OmegaBackground"]
            end
            cached_messages[nodeindex] = Dict()
            parent = x.parent
            x.parent = nothing
            for (alpha, omega) in alpha_and_single_omega_grid
                model_time += @elapsed model = Omega_model_func(cached_model,omega,alpha,GTRmat,F3x4_freqs,code)
                felsenstein!(x, model)
                copy_time += @elapsed cached_messages[nodeindex][(alpha, omega)] = deepcopy(x.message)
            end
            x.parent = parent
            x.children = FelNode[]
        end
    end

    @time @show model_time copy_time Base.summarysize(cached_messages) / (1024^2) "mb" length(pure_subclades) Base.summarysize(cached_model) / (1024^2) "mb"
    #return pure_subclades, cached_messages, cached_tag_inds
    for (row_ind,cp) in enumerate(codon_param_vec)
        alpha = cp[1]
        omegas = cp[2:end]
        tagged_models = N_Omegas_model_func(cached_model,tags,omegas,alpha,GTRmat,F3x4_freqs, code)

        for x in pure_subclades
            nodeindex = x.nodeindex
            x.message = cached_messages[nodeindex][(alpha, omegas[cached_tag_inds[nodeindex]])]
        end

        felsenstein!(tree,tagged_models)
        #This combine!() is needed because the current site_LLs function applies to a partition
        #And after a felsenstein pass, you don't have the eq freqs factored in.
        #We could make a version of log_likelihood() that returns the partitions instead of just the sum
        combine!.(tree.message,tree.parent_message)

        log_con_lik_matrix[row_ind,:] .= MolecularEvolution.site_LLs(tree.message[1]) #Check that these grab the scaling constants as well!
        verbosity > 0 && if mod(row_ind,500)==1
            print(round(100*row_ind/length(codon_param_vec)),"% ")
            flush(stdout)
        end
    end

    verbosity > 0 && println()

    con_lik_matrix = zeros(size(log_con_lik_matrix));
    site_scalers = maximum(log_con_lik_matrix, dims = 1);
    for i in 1:num_sites
        con_lik_matrix[:,i] .= exp.(log_con_lik_matrix[:,i] .- site_scalers[i])
    end

    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end

function difFUBAR_grid_tree_surgery_chunks(tree, tags, GTRmat, F3x4_freqs, code; verbosity = 1, foreground_grid = 6, background_grid = 4)

    MolecularEvolution.set_node_indices!(tree)
    @time trees = [tree, [deepcopy(tree) for _ = 1:(Threads.nthreads() - 1)]...]
    cached_models = [MG94_cacher(code) for _ = 1:Threads.nthreads()]
    @time nodelists = [getnodelist(tree) for tree in trees]

    #Defines the grid used for inference.
    function gridsetup(lb, ub, num_below_one, trin, tr)
        step = (trin(1.0) - trin(lb))/num_below_one
        return tr.(trin(lb):step:trin(ub))
    end
    tr(x) = 10^x-0.05
    trinv(x) =  log10(x+0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr); 
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree) if !isroot(n)]) > num_groups
    tensor_dims = 1+num_groups+is_background;
    
    codon_param_vec = [[a] for a in alphagrid]
    param_kinds = ["Alpha"]
    for g in 1:num_groups
        push!(param_kinds, "OmegaG$(g)")
        codon_param_vec = add_to_each_element(codon_param_vec,omegagrid)
    end
    if is_background
        push!(param_kinds, "OmegaBackground")
        codon_param_vec = add_to_each_element(codon_param_vec,background_omega_grid)
    end
    codon_param_vec;
    
    pure_subclades, _, _ = getpuresubclades(tree, tags)

    if length(pure_subclades) > 0
        alpha_and_single_omega_grids = Dict()
        alphagrid_vectorized = [[a] for a in alphagrid]
        alpha_and_single_omega_grids["Omega"] = add_to_each_element(alphagrid_vectorized,omegagrid)
        if is_background
            alpha_and_single_omega_grids["OmegaBackground"] = add_to_each_element(alphagrid_vectorized,background_omega_grid)
        end
    end

    function do_subgrid(pure_subclade, cached_model, aso_chunk)
        cached_messages_x = Dict()
        for (alpha, omega) in aso_chunk
            model = Omega_model_func(cached_model,omega,alpha,GTRmat,F3x4_freqs,code)
            felsenstein!(pure_subclade, model)
            cached_messages_x[(alpha, omega)] = deepcopy(pure_subclade.message)
        end
        return cached_messages_x
    end

    cached_messages = Dict()
    cached_tag_inds = Dict()
    for x in pure_subclades
        tag_ind_below = model_ind(x.children[1].name, tags)
        nodeindex = x.nodeindex
        cached_tag_inds[nodeindex] = tag_ind_below
        if tag_ind_below <= num_groups
            alpha_and_single_omega_grid = alpha_and_single_omega_grids["Omega"]
        else
            alpha_and_single_omega_grid = alpha_and_single_omega_grids["OmegaBackground"]
        end

        parents = FelNode[]
        for nodelist in nodelists
            push!(parents, nodelist[nodeindex].parent)
            nodelist[nodeindex].parent = nothing
        end

        aso_chunks = Iterators.partition(alpha_and_single_omega_grid, max(1, ceil(Int, length(alpha_and_single_omega_grid) / Threads.nthreads())))
        tasks = []
        for (i, aso_chunk) in enumerate(aso_chunks)
            # Spawn the task and add it to the array
            task = Threads.@spawn do_subgrid(nodelists[i][nodeindex], cached_models[i], aso_chunk)
            push!(tasks, task)
        end

        # Wait for all tasks to finish and collect their return values
        @time dicts = [fetch(task) for task in tasks];

        # Merge all the returned Dicts, and put it in the big Dict of messages
        merged_dicts = merge(dicts...)
        cached_messages[nodeindex] = merged_dicts

        for (nodelist, parent) in zip(nodelists, parents)
            nodelist[nodeindex].parent = parent
            nodelist[nodeindex].children = FelNode[]
        end
    end

    GC.gc()

    num_sites = tree.message[1].sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l,num_sites);

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")

    function do_subgrid(tree, cached_model, cpv_chunk, idx)
        # Note that cpv_chunk is already enumerated
        for (row_ind,cp) in cpv_chunk
            alpha = cp[1]
            omegas = cp[2:end]
            tagged_models = N_Omegas_model_func(cached_model, tags,omegas,alpha,GTRmat,F3x4_freqs, code)

            for x in pure_subclades
                nodeindex = x.nodeindex
                # Get the local equivalent node to x
                y = nodelists[idx][nodeindex]
                y.message = cached_messages[nodeindex][(alpha, omegas[cached_tag_inds[nodeindex]])]
            end

            felsenstein!(tree,tagged_models)
            #This combine!() is needed because the current site_LLs function applies to a partition
            #And after a felsenstein pass, you don't have the eq freqs factored in.
            #We could make a version of log_likelihood() that returns the partitions instead of just the sum
            combine!.(tree.message,tree.parent_message)
            log_con_lik_matrix[row_ind,:] .= MolecularEvolution.site_LLs(tree.message[1]) #Check that these grab the scaling constants as well!
            #verbosity > 0 && if mod(row_ind,500)==1
            #    print(round(100*row_ind/length(codon_param_vec)),"% ")
            #    flush(stdout)
            #end
        end
    end

    cpv_chunks = Iterators.partition(enumerate(codon_param_vec), max(1, ceil(Int, length(codon_param_vec) / Threads.nthreads())))
    tasks = []
    for (i, cpv_chunk) in enumerate(cpv_chunks)
        # Spawn the task and add it to the array
        task = Threads.@spawn do_subgrid(trees[i], cached_models[i], cpv_chunk, i)
        push!(tasks, task)
    end

    # Wait for all tasks to finish
    foreach(wait, tasks)

    verbosity > 0 && println()

    con_lik_matrix = zeros(size(log_con_lik_matrix));
    site_scalers = maximum(log_con_lik_matrix, dims = 1);
    for i in 1:num_sites
        con_lik_matrix[:,i] .= exp.(log_con_lik_matrix[:,i] .- site_scalers[i])
    end

    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end

function difFUBAR_sample(con_lik_matrix, iters; verbosity = 1)
    verbosity > 0 &&  println("Step 4: Running Gibbs sampler to infer site categories.")
    alloc_grid,theta = LDA_gibbs_track_allocation_vec(con_lik_matrix,0.1, iters = iters)
    return alloc_grid,theta
end

function difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity = 1, sites_to_plot = nothing, exports = true)
    grid_size, num_sites = size(alloc_grid)

    r(s) = round(s,digits = 4);

    detected_sites = Int64[]
    group1_volumes = Vector{Float64}[]
    group2_volumes = Vector{Float64}[]
    alpha_volumes = Vector{Float64}[]
    detections = Vector{Float64}[] #legacy name - now includes all 4 "relevant" site posteriors
    param_means = Vector{Float64}[]

    ω1 = [c[2] for c in codon_param_vec];
    ω2 = [c[3] for c in codon_param_vec];
    alphas = [c[1] for c in codon_param_vec];
    ω1_greater_filt = ω1 .> ω2;
    ω2_greater_filt = ω2 .> ω1;
    ω1_pos_filt = ω1 .> 1.0;
    ω2_pos_filt = ω2 .> 1.0;

    verbosity > 0 && println("Step 5: Tabulating and plotting. Detected sites:")
    for site in 1:num_sites
        ω1_greater_posterior = sum(alloc_grid[ω1_greater_filt,site])/sum(alloc_grid[:,site])
        ω2_greater_posterior = sum(alloc_grid[ω2_greater_filt,site])/sum(alloc_grid[:,site])
        ω1_pos_posterior = sum(alloc_grid[ω1_pos_filt,site])/sum(alloc_grid[:,site])
        ω2_pos_posterior = sum(alloc_grid[ω2_pos_filt,site])/sum(alloc_grid[:,site])
        detecs = [ω1_greater_posterior,ω2_greater_posterior,ω1_pos_posterior,ω2_pos_posterior]
        
        site_counts_ω1 = collapse_counts(ω1,alloc_grid[:,site], cases = omegagrid)
        site_counts_ω2 = collapse_counts(ω2,alloc_grid[:,site], cases = omegagrid)
        site_counts_alphas = collapse_counts(alphas,alloc_grid[:,site], cases = alphagrid)
        
        mean_alpha = sum(site_counts_alphas .* alphagrid)
        mean_ω1 = sum(site_counts_ω1 .* omegagrid)
        mean_ω2 = sum(site_counts_ω2 .* omegagrid)
        
        push!(detections,detecs)
        push!(param_means,[mean_alpha,mean_ω1,mean_ω2])
        push!(group1_volumes,site_counts_ω1)
        push!(group2_volumes,site_counts_ω2)
        push!(alpha_volumes,site_counts_alphas)
        
        if maximum(detecs) > pos_thresh
            verbosity > 0 && print("Site $(site) - ")
            verbosity > 0 && print("P(ω1 > ω2):", ω1_greater_posterior)
            verbosity > 0 && print("; P(ω2 > ω1):", ω2_greater_posterior)
            verbosity > 0 && print("; P(ω1 > 1):", ω1_pos_posterior)
            verbosity > 0 && println("; P(ω2 > 1):", ω2_pos_posterior)
            push!(detected_sites,site)
        end
    end

    #Exporting site data
    df = DataFrame()
    df[!,"Codon Sites"] = [1:num_sites;]
    df[!,"P(ω1 > ω2)"] = [d[1] for d in detections]
    df[!,"P(ω2 > ω1)"] = [d[2] for d in detections]
    df[!,"P(ω1 > 1)"] = [d[3] for d in detections]
    df[!,"P(ω2 > 1)"] = [d[4] for d in detections];
    df[!,"mean(α)"] = [d[1] for d in param_means]
    df[!,"mean(ω1)"] = [d[2] for d in param_means]
    df[!,"mean(ω2)"] = [d[3] for d in param_means]

    verbosity > 0 && println("\nIf exports = true, writing results for all sites to CSV: "*analysis_name*"_posteriors.csv")
    exports && CSV.write(analysis_name*"_posteriors.csv", df)

    sites = [1:num_sites;]

    #Select the sites that will get plotted, in case you want to customize this.
    if isnothing(sites_to_plot)
        sites_to_plot = detected_sites
    end

    if length(sites_to_plot) == 0
        verbosity > 0 && println("No sites detected above threshold.")
    elseif exports
        verbosity > 0 && println("Plotting alpha and omega distributions. If exports = true, saved as "*analysis_name*"_violin_*.pdf")
    
        #Assumes alpha and omega grids are the same!? Currently enforced by args passed into difFUBAR_grid
        #Maybe this is ok
        grd = round.(omegagrid, digits = 3) 

        #Three plotting examples.
        #Plot the alphas for each flagged site

        figure(figsize = (3,1+length(sites[sites_to_plot])/2))
        FUBAR_violin_plot(sites[sites_to_plot],alpha_volumes[sites_to_plot] .* 0.75,grd, tag = "α", color = "green", x_label = "α")
        savefig(analysis_name*"_violin_alpha.pdf");
        close()

        #Plot the G1 and G2 omegas
        figure(figsize = (3,1+length(sites[sites_to_plot])/2))
        FUBAR_violin_plot(sites[sites_to_plot],group1_volumes[sites_to_plot],grd, tag = "ω1", color = tag_colors[1])
        FUBAR_violin_plot(sites[sites_to_plot],group2_volumes[sites_to_plot],grd, tag = "ω2", color = tag_colors[2], x_label = "ω")
        savefig(analysis_name*"_violin_omegas.pdf");
        close()

        #Plot all three parameters, using the v_offset to separate the alphas from the omegas
        figure(figsize = (3,4+length(sites[sites_to_plot])/2))
        FUBAR_violin_plot(sites[sites_to_plot],group1_volumes[sites_to_plot] .* 0.5,grd, tag = "ω1", color = tag_colors[1], v_offset = -0.1)
        FUBAR_violin_plot(sites[sites_to_plot],group2_volumes[sites_to_plot] .* 0.5,grd, tag = "ω2", color = tag_colors[2], v_offset = -0.1)
        FUBAR_violin_plot(sites[sites_to_plot],alpha_volumes[sites_to_plot] .* 0.5,grd, tag = "α", color = "green", v_offset = 0.1)
        savefig(analysis_name*"_violin_all_params.pdf");
        close()

        #Coerce the violin plot function to also viz the "detection" posteriors.
        figure(figsize = (1.5,1+length(sites[sites_to_plot])/2))
        floored_detec = [clamp.((d .- 0.95).*20 , 0.0, 1.0) for d in detections[sites_to_plot]]
        FUBAR_violin_plot(sites[sites_to_plot],[[f[1],0.0,0.0,0.0] for f in floored_detec].*0.5,
            ["P(ω1>ω2)","P(ω2>ω1)","P(ω1>1)","P(ω2>1)"], tag = "P(ω1>ω2)", color = tag_colors[1], 
            vertical_ind = nothing, plot_legend = false)
        FUBAR_violin_plot(sites[sites_to_plot],[[0.0,f[2],0.0,0.0] for f in floored_detec].*0.5,
            ["P(ω1>ω2)","P(ω2>ω1)","P(ω1>1)","P(ω2>1)"], tag = "P(ω2>ω1)", color = tag_colors[2], 
            vertical_ind = nothing, plot_legend = false)
        FUBAR_violin_plot(sites[sites_to_plot],[[0.0,0.0,f[3],0.0] for f in floored_detec].*0.5,
            ["P(ω1>ω2)","P(ω2>ω1)","P(ω1>1)","P(ω2>1)"], tag = "P(ω1>1)", color = tag_colors[1], 
            vertical_ind = nothing, plot_legend = false)
        FUBAR_violin_plot(sites[sites_to_plot],[[0.0,0.0,0.0,f[4]] for f in floored_detec].*0.5,
            ["P(ω1>ω2)","P(ω2>ω1)","P(ω1>1)","P(ω2>1)"], tag = "P(ω2>1)", color = tag_colors[2], 
            vertical_ind = nothing, legend_ncol = 2, x_label = "", plot_legend = false)
        savefig(analysis_name*"_detections.pdf");
        close()

        
    end
    
    if exports
        #A plot of the omega means for all sites.
        figure(figsize = (20,3))
        omega1_means = [p[2] for p in param_means] 
        omega2_means = [p[3] for p in param_means]
        for i in 1:length(omega1_means)
            tc = "black"
            if omega1_means[i] > omega2_means[i]
                tc = tag_colors[1]
            else
                tc = tag_colors[2]
            end
            
            diff_mul = 1.0
            if !(maximum(detections[i][1:2]) > pos_thresh)
                diff_mul = 0.15
            end
            
            pos1_mul = 1.0
            if !(detections[i][3] > pos_thresh)
                pos1_mul = 0.2
            end
            
            pos2_mul = 1.0
            if !(detections[i][4] > pos_thresh)
                pos2_mul = 0.2
            end

            plot([i,i],[omega1_means[i],omega2_means[i]], color = tc, alpha = 0.75*diff_mul)
            plot([i], [omega1_means[i]],".", color = tag_colors[1], alpha = 0.75*pos1_mul)
            plot([i], [omega2_means[i]],".", color = tag_colors[2], alpha = 0.75*pos2_mul)
        end

        plot([-100,-100],[2,2],".",color = tag_colors[1], alpha = 0.75, label = "ω1>1")
        plot([-100,-100],[2,2],".",color = tag_colors[2], alpha = 0.75, label = "ω2>1")
        plot([-100,-100],[2,2], color = tag_colors[1], alpha = 0.75, label = "ω1>ω2")
        plot([-100,-100],[2,2], color = tag_colors[2], alpha = 0.75, label = "ω2>ω1")
        plot([-2,num_sites+3], [1.0,1.0],"--", color = "grey", alpha = 0.5, label = "ω=1")
        xlabel("Codon Sites")
        ylabel("ω")
        yscale("symlog")
        yticks([0.01, 0.5, 1.0, 1.6, 2.5, 5.0, 10.0],[0.01, 0.5, 1.0, 1.6, 2.5, 5.0, 10.0])
        xlim(-4,num_sites+5)
        legend(loc="upper center", bbox_to_anchor=(0.5, 1.15), ncol=5, fancybox=true, shadow=true)
        tight_layout()
        savefig(analysis_name*"_site_omega_means.pdf");
    end

    return df
end
export difFUBAR_tabulate

#Must return enough to re-calculate detections etc
function difFUBAR(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh = 0.95, iters = 2500, verbosity = 1, exports = true, code = MolecularEvolution.universal_code)
    analysis_name = outpath
    tree, tags, tag_colors, analysis_name = difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports = exports, verbosity = verbosity)
    tree, alpha,beta,GTRmat,F3x4_freqs,eq_freqs = difFUBAR_global_fit(seqnames, seqs, tree, generate_tag_stripper(tags), code, verbosity = verbosity)
    con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, 
                                                                                verbosity = verbosity, foreground_grid = 6, background_grid = 4)
    alloc_grid,theta = difFUBAR_sample(con_lik_matrix, iters, verbosity = verbosity)
    df = difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity = verbosity, exports = exports)

    #Return df, (tuple of partial calculations needed to re-run tablulate)
    return df, (alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors)
end
export difFUBAR



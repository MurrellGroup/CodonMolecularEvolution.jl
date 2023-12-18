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

# added for parallellizing difFUBAR_grid

function model_ind(str::String, tags::Vector{String})
    ind = length(tags) + 1
    for (i, t) in enumerate(tags)
        if occursin(t, str)
            ind = i
        end
    end
    return ind
end

function collapse_counts(param_vec, count_vec; cases=nothing)
    if isnothing(cases)
        cases = sort(union(param_vec))
    end
    d = Dict(zip(cases, 1:length(cases)))
    storage = zeros(Int64, length(cases))
    for i in 1:length(count_vec)
        storage[d[param_vec[i]]] += count_vec[i]
    end
    return storage ./ sum(storage)
end

function FUBAR_violin_plot(sites, group1_volumes, omegagrid;
    color="black", tag="", alpha=0.6,
    x_label="Parameter", y_label="Codon Sites",
    v_offset=0.0, legend_ncol=3,
    vertical_ind=findfirst(omegagrid .>= 1.0),
    plot_legend=true)

    ypos = [-i * 0.5 for i in 1:length(sites)]
    if !isnothing(vertical_ind)
        bar!([vertical_ind], [2 + maximum(ypos) - minimum(ypos)], bottom=[minimum(ypos) - 1], color="grey", alpha=0.05, label=:none)
    end

    for i in 1:length(sites)
        center_line = ypos[i]
        a = 1:length(omegagrid)
        b = group1_volumes[i]
        c = (v_offset .+ center_line .+ (-0.5 .* group1_volumes[i]))

        bar!(a, b + c, fillto=c, linewidth=0, bar_edges=false, linealpha=0.0, ylims=(minimum(c) - 1, 0), color=color, alpha=alpha, label=:none)

    end

    bar!([1], [1], bottom=[1000], color=color, alpha=alpha, label=tag, linewidth=0, bar_edges=false, linealpha=0.0)
    bar!(yticks=(ypos, sites))
    bar!(xticks=((1:length(omegagrid)), omegagrid), xrotation=90)
    bar!(ylabel=y_label, xlabel=x_label)

    bar!(ylim=(minimum(ypos) - 0.5, maximum(ypos) + 0.5))
    if plot_legend
        plot!(
            legend=:outertop,
            legendcolumns=legend_ncol,
            shadow=true, fancybox=true, bbox_to_anchor=(0.5, 1.01)
        )
    end
end

function FUBAR_omega_plot(param_means, tag_colors, pos_thresh, detections, num_sites)
    #A plot of the omega means for all sites.
    omega1_means = [p[2] for p in param_means]
    omega2_means = [p[3] for p in param_means]

    t(x) = log10(x + 1)
    invt(y) = 10^y - 1

    omega1_means = t.(omega1_means)
    omega2_means = t.(omega2_means)

    for i in 1:length(omega1_means)
        tc = "black"
        if omega1_means[i] > omega2_means[i]
            tc = tag_colors[1]
        else
            tc = tag_colors[2]
        end

        diff_mul = 1.0
        if !(maximum(detections[i][1:2]) > pos_thresh)
            diff_mul = 0.1
        end

        pos1_mul = 1.0
        if !(detections[i][3] > pos_thresh)
            pos1_mul = 0.15
        end

        pos2_mul = 1.0
        if !(detections[i][4] > pos_thresh)
            pos2_mul = 0.15
        end
        plot!([i, i], [omega1_means[i], omega2_means[i]], color=tc, alpha=0.75 * diff_mul, linewidth=2, xlim=(-4, num_sites + 5), label="", yscale=:log10)
        scatter!([i], [omega1_means[i]], color=tag_colors[1], alpha=0.75 * pos1_mul, ms=2.5, label="", markerstrokecolor=:auto, yscale=:log10)
        scatter!([i], [omega2_means[i]], color=tag_colors[2], alpha=0.75 * pos2_mul, ms=2.5, label="", markerstrokecolor=:auto, yscale=:log10)
    end


    scatter!([-100, -100], [2, 2], color=tag_colors[1], alpha=0.75, label="ω1>1", ms=2.5)
    scatter!([-100, -100], [2, 2], color=tag_colors[2], alpha=0.75, label="ω2>1", ms=2.5)
    plot!([-100, -100], [2, 2], color=tag_colors[1], alpha=0.75, label="ω1>ω2", linewidth=2)
    plot!([-100, -100], [2, 2], color=tag_colors[2], alpha=0.75, label="ω2>ω1", linewidth=2)
    plot!([-2, num_sites + 3], [log10(1.0 + 1), log10(1.0 + 1)], color="grey", alpha=0.5, label="ω=1", linestyle=:dot, linewidth=2)
    xlabel!("Codon Sites")
    ylabel!("ω")

    n_points = 8
    lb = 0.01
    ub = 10
    points = collect(t(lb):(t(ub)-t(lb))/(n_points-1):t(ub))
    ticklabels = string.(round.(invt.(points), sigdigits=2))
    yticks!(points, ticklabels)

    xticks!(0:50:num_sites)

    plot!(
        legend=:outertop,
        legendcolumns=5,
        ylim=(0, log10(11)))

end

"""


"""
function difFUBAR_init(outpath_and_file_prefix, treestring, tags, tag_colors; verbosity=1, exports=true, strip_tags_from_name=generate_tag_stripper(tags))

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
    verbosity > 0 && println("Step 1: Initialization. If exports = true, tree showing the assignment of branches to groups/colors will be exported to: " * analysis_name * "_tagged_input_tree.svg.")


    p = sortperm(tags)
    tags, tag_colors = tags[p], tag_colors[p]
    push!(tag_colors, "black") #ASSUMPTION: background color is always black

    if exports
        #Replace with Phylo.jl based plot?
        color_dict = Dict(zip(getnodelist(tree), [tag_colors[model_ind(n.name, tags)] for n in getnodelist(tree)]))
        label_dict = Dict(zip(getnodelist(tree), [strip_tags_from_name(n.name) for n in getnodelist(tree)]))
        img = tree_draw(tree, canvas_height=(3 + length(getleaflist(tree)) / 5)cm,
            draw_labels=true, dot_color_dict=color_dict,
            line_color_dict=color_dict, line_width=0.3, min_dot_size=0.01,
            nodelabel_dict=label_dict)
        img |> SVG(analysis_name * "_tagged_input_tree.svg", 15cm, (3 + length(getleaflist(tree)) / 5)cm)
    end

    #Tags and tag colors are now ordered, and tag_colors includes the untagged category
    return tree, tags, tag_colors, analysis_name
end

function difFUBAR_global_fit(seqnames, seqs, tree, leaf_name_transform, code; verbosity=1)

    verbosity > 0 && println("Step 2: Optimizing global codon model parameters.")

    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = optimize_MG94_F3x4(seqnames, seqs, tree, leaf_name_transform=leaf_name_transform, genetic_code=code)

    ######
    #optionally polish branch lengths and topology
    ######

    return tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs
end

#foreground_grid and background_grid control the number of categories below 1.0
function difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code; verbosity=1, foreground_grid=6, background_grid=4)

    cached_model = MG94_cacher(code)

    #This is the function that assigns models to branches
    #Sometimes there will be the same number of tags as omegas, but sometimes there will be one more omega.
    #Depending on whether there are any background branches (UNTESTED WITH NO BACKGROUND BRANCHES)
    function N_Omegas_model_func(tags, omega_vec, alpha, nuc_mat, F3x4, code)
        models = [cached_model(alpha, alpha * o, nuc_mat, F3x4, genetic_code=code) for o in omega_vec]
        return n::FelNode -> [models[model_ind(n.name, tags)]]
    end

    #Defines the grid used for inference.
    function gridsetup(lb, ub, num_below_one, trin, tr)
        step = (trin(1.0) - trin(lb)) / num_below_one
        return tr.(trin(lb):step:trin(ub))
    end
    tr(x) = 10^x - 0.05
    trinv(x) = log10(x + 0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree)]) > num_groups
    #is_background = false
    tensor_dims = 1 + num_groups + is_background

    function add_to_each_element(vec_of_vec, elems)
        return [vcat(v, [e]) for v in vec_of_vec for e in elems]
    end

    codon_param_vec = [[a] for a in alphagrid]
    param_kinds = ["Alpha"]
    for g in 1:num_groups
        push!(param_kinds, "OmegaG$(g)")
        codon_param_vec = add_to_each_element(codon_param_vec, omegagrid)
    end
    if is_background
        push!(param_kinds, "OmegaBackground")
        codon_param_vec = add_to_each_element(codon_param_vec, background_omega_grid)
    end
    codon_param_vec

    num_sites = tree.message[1].sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l, num_sites)

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")


    # modify the tree 

    # compute recursion G1
    # compute recursion G2

    # How does the omega influence the propagation?


    for (row_ind, cp) in enumerate(codon_param_vec)
        alpha = cp[1]
        omegas = cp[2:end]
        tagged_models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)

        #print("felsenstein")
        felsenstein!(tree, tagged_models)
        #This combine!() is needed because the current site_LLs function applies to a partition
        #And after a felsenstein pass, you don't have the eq freqs factored in.
        #We could make a version of log_likelihood() that returns the partitions instead of just the sum
        #print("combine")
        combine!.(tree.message, tree.parent_message)
        #println(length(tree.message[1]))
        #print("log_con_like_matrix")
        log_con_lik_matrix[row_ind, :] .= MolecularEvolution.site_LLs(tree.message[1]) #Check that these grab the scaling constants as well!
        verbosity > 0 && if mod(row_ind, 500) == 1
            print(round(100 * row_ind / length(codon_param_vec)), "% ")
            flush(stdout)
        end
    end

    verbosity > 0 && println()

    con_lik_matrix = zeros(size(log_con_lik_matrix))
    site_scalers = maximum(log_con_lik_matrix, dims=1)
    for i in 1:num_sites
        con_lik_matrix[:, i] .= exp.(log_con_lik_matrix[:, i] .- site_scalers[i])
    end

    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end

#foreground_grid and background_grid control the number of categories below 1.0
function difFUBAR_grid_pruned_1(tree, tags, GTRmat, F3x4_freqs, code; verbosity=1, foreground_grid=6, background_grid=4)

    cached_model = MG94_cacher(code)

    #This is the function that assigns models to branches
    #Sometimes there will be the same number of tags as omegas, but sometimes there will be one more omega.
    #Depending on whether there are any background branches (UNTESTED WITH NO BACKGROUND BRANCHES)
    function N_Omegas_model_func(tags, omega_vec, alpha, nuc_mat, F3x4, code)
        models = [cached_model(alpha, alpha * o, nuc_mat, F3x4, genetic_code=code) for o in omega_vec]
        return n::FelNode -> [models[model_ind(n.name, tags)]]
    end

    #Defines the grid used for inference.
    function gridsetup(lb, ub, num_below_one, trin, tr)
        step = (trin(1.0) - trin(lb)) / num_below_one
        return tr.(trin(lb):step:trin(ub))
    end
    tr(x) = 10^x - 0.05
    trinv(x) = log10(x + 0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    #is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree)]) > num_groups
    is_background = false
    tensor_dims = 1 + num_groups + is_background

    function add_to_each_element(vec_of_vec, elems)
        return [vcat(v, [e]) for v in vec_of_vec for e in elems]
    end

    codon_param_vec = [[a] for a in alphagrid]
    param_kinds = ["Alpha"]
    for g in 1:num_groups
        push!(param_kinds, "OmegaG$(g)")
        codon_param_vec = add_to_each_element(codon_param_vec, omegagrid)
    end
    if is_background
        push!(param_kinds, "OmegaBackground")
        codon_param_vec = add_to_each_element(codon_param_vec, background_omega_grid)
    end
    codon_param_vec

    num_sites = tree.message[1].sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l, num_sites)

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")


    # modify the tree 

    # compute recursion G1
    # compute recursion G2

    # How does the omega influence the propagation?

    pure_subclades = FelNode[]
    # Traverses the tree recursively with a dfs whilst pushing roots of pure subclades to list^
    function find_pure_subclades(node)
        # Get the index of the node's tag
        tag_ind_of_node = model_ind(node.name, tags)

        # If the node is a leaf, it's pure
        if length(node.children) == 0
            return true, tag_ind_of_node
        end

        children_is_pure_and_tag = []
        pure_children = FelNode[]

        for child in node.children
            child_is_pure, tag_ind = find_pure_subclades(child)
            if child_is_pure
                push!(pure_children, child)
            end
            push!(children_is_pure_and_tag, (child_is_pure, tag_ind))
        end

        # Get the index of the node's first child's tag
        tag_ind_of_first_child = children_is_pure_and_tag[1][2]

        # This is the case where the subclade starting at node is pure
        if all(x == (true, tag_ind_of_first_child) for x in children_is_pure_and_tag)
            if tag_ind_of_node != tag_ind_of_first_child
                # The purity is broken at this node
                push!(pure_subclades, node)
                return false, tag_ind_of_node
            end
            # The purity is not broken at this node
            return true, tag_ind_of_node
        end

        # This is the case where some child has mixed tags or the children are pure with regards to different tags
        for pure_child in pure_children
            if length(pure_child.children) == 0
                # We don't want to push leaves into pure_subclades
                continue
            end
            push!(pure_subclades, pure_child)
        end
        return false, tag_ind_of_node
    end

    time_elapsed = @elapsed find_pure_subclades(tree)
    println("Max prune " * string(time_elapsed))

    cached_messages = Dict()
    cached_tag_inds = Dict()
    for x in pure_subclades
        cached_messages[x] = Dict()
        cached_tag_inds[x] = model_ind(x.children[1].name, tags)
        parent = x.parent
        x.parent = nothing
        for cp in codon_param_vec
            alpha = cp[1]
            omegas = cp[2:end]

            relevant_omega = omegas[cached_tag_inds[x]]

            if haskey(cached_messages[x], (alpha, relevant_omega))
                continue
            end

            models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)
            felsenstein!(x, models)
            cached_messages[x][(alpha, relevant_omega)] = deepcopy(x.message)
        end
        x.parent = parent
        x.children = FelNode[]
    end

    for (row_ind, cp) in enumerate(codon_param_vec)
        alpha = cp[1]
        omegas = cp[2:end]
        tagged_models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)

        # added for subtree prune 
        for x in pure_subclades
            x.message = cached_messages[x][(alpha, omegas[cached_tag_inds[x]])]
        end

        felsenstein!(tree, tagged_models)
        #This combine!() is needed because the current site_LLs function applies to a partition
        #And after a felsenstein pass, you don't have the eq freqs factored in.
        #We could make a version of log_likelihood() that returns the partitions instead of just the sum
        combine!.(tree.message, tree.parent_message)

        log_con_lik_matrix[row_ind, :] .= MolecularEvolution.site_LLs(tree.message[1]) #Check that these grab the scaling constants as well!
        verbosity > 0 && if mod(row_ind, 500) == 1
            print(round(100 * row_ind / length(codon_param_vec)), "% ")
            flush(stdout)
        end
    end

    verbosity > 0 && println()

    con_lik_matrix = zeros(size(log_con_lik_matrix))
    site_scalers = maximum(log_con_lik_matrix, dims=1)
    for i in 1:num_sites
        con_lik_matrix[:, i] .= exp.(log_con_lik_matrix[:, i] .- site_scalers[i])
    end

    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end



#foreground_grid and background_grid control the number of categories below 1.0
function difFUBAR_grid_pruned_2(tree, tags, GTRmat, F3x4_freqs, code; verbosity=1, foreground_grid=6, background_grid=4)

    cached_model = MG94_cacher(code)

    #This is the function that assigns models to branches
    #Sometimes there will be the same number of tags as omegas, but sometimes there will be one more omega.
    #Depending on whether there are any background branches (UNTESTED WITH NO BACKGROUND BRANCHES)
    function N_Omegas_model_func(tags, omega_vec, alpha, nuc_mat, F3x4, code)
        models = [cached_model(alpha, alpha * o, nuc_mat, F3x4, genetic_code=code) for o in omega_vec]
        return n::FelNode -> [models[model_ind(n.name, tags)]]
    end

    #Defines the grid used for inference.
    function gridsetup(lb, ub, num_below_one, trin, tr)
        step = (trin(1.0) - trin(lb)) / num_below_one
        return tr.(trin(lb):step:trin(ub))
    end
    tr(x) = 10^x - 0.05
    trinv(x) = log10(x + 0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree)]) > num_groups
    tensor_dims = 1 + num_groups + is_background

    function add_to_each_element(vec_of_vec, elems)
        return [vcat(v, [e]) for v in vec_of_vec for e in elems]
    end

    codon_param_vec = [[a] for a in alphagrid]
    param_kinds = ["Alpha"]
    for g in 1:num_groups
        push!(param_kinds, "OmegaG$(g)")
        codon_param_vec = add_to_each_element(codon_param_vec, omegagrid)
    end
    if is_background
        push!(param_kinds, "OmegaBackground")
        codon_param_vec = add_to_each_element(codon_param_vec, background_omega_grid)
    end
    codon_param_vec

    num_sites = tree.message[1].sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l, num_sites)

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")


    # modify the tree 

    # compute recursion G1
    # compute recursion G2

    # How does the omega influence the propagation?

    # getnodelist and get a specific node from nodelist
    #function get_node_by_name(query_node::String, tree::FelNode)
    #    # This is a useful function that maybe should be in AbstractTreeNode.jl
    #    for node in getnodelist(tree)
    #        if node.name == query_node
    #            return node
    #        end
    #    end
    #end

    function check_purity_from_node_and_forward_in_time(tree)
        # collect all nodes under and see if they are the same group
        # if return true that means everything under the tree is pure (the tree node can be differnt though)
        node_groups = []
        for node in getnodelist(tree)
            if tree.name != node.name # we only check nodes under the input tree or subtree
                push!(node_groups, model_ind(node.name, tags))
            end
        end
        if length(unique(node_groups)) == 1
            return true
        else
            return false
        end
    end

    pure_clades = []
    function traverse_tree_to_check_for_pure_clades(pure_clades, tree)
        for child in tree.children
            if check_purity_from_node_and_forward_in_time(child)
                #print(child.name)
                #println(" PURE")
                # add this node, as everything under this node is pure
                push!(pure_clades, child)
            else
                #print(child.name)
                #println(" not pure") # keep recursion to se if we go forward in time if we can find pure clades
                traverse_tree_to_check_for_pure_clades(pure_clades, child)
            end
        end
        return pure_clades
    end

    time_elapsed = @elapsed subtree_tops = traverse_tree_to_check_for_pure_clades(pure_clades, tree)
    println("Patrick prune " * string(time_elapsed))

    cached_values = Dict()
    for x in subtree_tops
        #x = get_node_by_name(x, tree)
        cached_values[x] = Dict()
        parent = x.parent
        x.parent = nothing
        for cp in codon_param_vec
            alpha = cp[1]
            omegas = cp[2:end]

            relevant_omega = omegas[model_ind(x.name, tags)] #here we need to change to child omega

            if haskey(cached_values[x], (alpha, relevant_omega))
                continue
            end

            models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)
            felsenstein!(x, models)
            cached_values[x][(alpha, relevant_omega)] = deepcopy(x.message)
        end
        x.parent = parent
        x.children = FelNode[]
    end
    #println("DEBUG")
    #println(pure_clades)

    for (row_ind, cp) in enumerate(codon_param_vec)
        alpha = cp[1]
        omegas = cp[2:end]
        tagged_models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)

        # added for subtree prune 
        for x in pure_clades
            #x = get_node_by_name(x, tree)
            x.message = cached_values[x][(alpha, omegas[model_ind(x.name, tags)])]
        end

        felsenstein!(tree, tagged_models)
        #This combine!() is needed because the current site_LLs function applies to a partition
        #And after a felsenstein pass, you don't have the eq freqs factored in.
        #We could make a version of log_likelihood() that returns the partitions instead of just the sum
        combine!.(tree.message, tree.parent_message)

        log_con_lik_matrix[row_ind, :] .= MolecularEvolution.site_LLs(tree.message[1]) #Check that these grab the scaling constants as well!
        verbosity > 0 && if mod(row_ind, 500) == 1
            print(round(100 * row_ind / length(codon_param_vec)), "% ")
            flush(stdout)
        end
    end

    verbosity > 0 && println()

    con_lik_matrix = zeros(size(log_con_lik_matrix))
    site_scalers = maximum(log_con_lik_matrix, dims=1)
    for i in 1:num_sites
        con_lik_matrix[:, i] .= exp.(log_con_lik_matrix[:, i] .- site_scalers[i])
    end

    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end

function difFUBAR_grid_pruned_3(tree, tags, GTRmat, F3x4_freqs, code; verbosity=1, foreground_grid=6, background_grid=4)

    cached_model = MG94_cacher(code)

    #This is the function that assigns models to branches
    #Sometimes there will be the same number of tags as omegas, but sometimes there will be one more omega.
    #Depending on whether there are any background branches (UNTESTED WITH NO BACKGROUND BRANCHES)
    function N_Omegas_model_func(tags, omega_vec, alpha, nuc_mat, F3x4, code)
        models = [cached_model(alpha, alpha * o, nuc_mat, F3x4, genetic_code=code) for o in omega_vec]
        return n::FelNode -> [models[model_ind(n.name, tags)]]
    end

    #Defines the grid used for inference.
    function gridsetup(lb, ub, num_below_one, trin, tr)
        step = (trin(1.0) - trin(lb)) / num_below_one
        return tr.(trin(lb):step:trin(ub))
    end
    tr(x) = 10^x - 0.05
    trinv(x) = log10(x + 0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    #is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree)]) > num_groups
    is_background = false
    tensor_dims = 1 + num_groups + is_background

    function add_to_each_element(vec_of_vec, elems)
        return [vcat(v, [e]) for v in vec_of_vec for e in elems]
    end

    codon_param_vec = [[a] for a in alphagrid]
    param_kinds = ["Alpha"]
    for g in 1:num_groups
        push!(param_kinds, "OmegaG$(g)")
        codon_param_vec = add_to_each_element(codon_param_vec, omegagrid)
    end
    if is_background
        push!(param_kinds, "OmegaBackground")
        codon_param_vec = add_to_each_element(codon_param_vec, background_omega_grid)
    end
    codon_param_vec

    num_sites = tree.message[1].sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l, num_sites)

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")


    # modify the tree 

    # compute recursion G1
    # compute recursion G2

    # How does the omega influence the propagation?

    # getnodelist and get a specific node from nodelist
    #function get_node_by_name(query_node::String, tree::FelNode)
    #    # This is a useful function that maybe should be in AbstractTreeNode.jl
    #    for node in getnodelist(tree)
    #        if node.name == query_node
    #            return node
    #        end
    #    end
    #end

    function check_purity_from_node_and_forward_in_time(tree)
        # collect all nodes under and see if they are the same group
        # if return true that means everything under the tree is pure (the tree node can be differnt though)
        node_groups = []
        for node in getnodelist(tree)
            if tree.name != node.name # we only check nodes under the input tree or subtree
                push!(node_groups, model_ind(node.name, tags))
            end
        end
        if length(unique(node_groups)) == 1
            return true
        else
            return false
        end
    end

    pure_clades = []
    function traverse_tree_to_check_for_pure_clades(pure_clades, tree)
        for child in tree.children
            if check_purity_from_node_and_forward_in_time(child)
                #print(child.name)
                #println(" PURE")
                # add this node, as everything under this node is pure
                push!(pure_clades, child)
            else
                #print(child.name)
                #println(" not pure") # keep recursion to se if we go forward in time if we can find pure clades
                traverse_tree_to_check_for_pure_clades(pure_clades, child)
            end
        end
        return pure_clades
    end

    time_elapsed = @elapsed subtree_tops = traverse_tree_to_check_for_pure_clades(pure_clades, tree)
    println("Patrick prune " * string(time_elapsed))

    cached_messages = Dict()
    cached_tag_inds = Dict()
    for x in subtree_tops
        cached_messages[x] = Dict()
        cached_tag_inds[x] = model_ind(x.children[1].name, tags)
        parent = x.parent
        x.parent = nothing
        for cp in codon_param_vec
            alpha = cp[1]
            omegas = cp[2:end]

            relevant_omega = omegas[cached_tag_inds[x]]

            if haskey(cached_messages[x], (alpha, relevant_omega))
                continue
            end

            models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)
            felsenstein!(x, models)
            cached_messages[x][(alpha, relevant_omega)] = deepcopy(x.message)
        end
        x.parent = parent
        x.children = FelNode[]
    end

    for (row_ind, cp) in enumerate(codon_param_vec)
        alpha = cp[1]
        omegas = cp[2:end]
        tagged_models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)

        # added for subtree prune 
        for x in subtree_tops
            x.message = cached_messages[x][(alpha, omegas[cached_tag_inds[x]])]
        end

        felsenstein!(tree, tagged_models)
        #This combine!() is needed because the current site_LLs function applies to a partition
        #And after a felsenstein pass, you don't have the eq freqs factored in.
        #We could make a version of log_likelihood() that returns the partitions instead of just the sum
        combine!.(tree.message, tree.parent_message)

        log_con_lik_matrix[row_ind, :] .= MolecularEvolution.site_LLs(tree.message[1]) #Check that these grab the scaling constants as well!
        verbosity > 0 && if mod(row_ind, 500) == 1
            print(round(100 * row_ind / length(codon_param_vec)), "% ")
            flush(stdout)
        end
    end

    verbosity > 0 && println()

    con_lik_matrix = zeros(size(log_con_lik_matrix))
    site_scalers = maximum(log_con_lik_matrix, dims=1)
    for i in 1:num_sites
        con_lik_matrix[:, i] .= exp.(log_con_lik_matrix[:, i] .- site_scalers[i])
    end

    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end


#foreground_grid and background_grid control the number of categories below 1.0
function difFUBAR_grid_pruned_4(tree, tags, GTRmat, F3x4_freqs, code; verbosity=1, foreground_grid=6, background_grid=4)

    cached_model = MG94_cacher(code)

    #This is the function that assigns models to branches
    #Sometimes there will be the same number of tags as omegas, but sometimes there will be one more omega.
    #Depending on whether there are any background branches (UNTESTED WITH NO BACKGROUND BRANCHES)
    function N_Omegas_model_func(tags, omega_vec, alpha, nuc_mat, F3x4, code)
        models = [cached_model(alpha, alpha * o, nuc_mat, F3x4, genetic_code=code) for o in omega_vec]
        return n::FelNode -> [models[model_ind(n.name, tags)]]
    end

    #Defines the grid used for inference.
    function gridsetup(lb, ub, num_below_one, trin, tr)
        step = (trin(1.0) - trin(lb)) / num_below_one
        return tr.(trin(lb):step:trin(ub))
    end
    tr(x) = 10^x - 0.05
    trinv(x) = log10(x + 0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    #is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree)]) > num_groups
    is_background = false

    tensor_dims = 1 + num_groups + is_background

    function add_to_each_element(vec_of_vec, elems)
        return [vcat(v, [e]) for v in vec_of_vec for e in elems]
    end

    codon_param_vec = [[a] for a in alphagrid]
    param_kinds = ["Alpha"]
    for g in 1:num_groups
        push!(param_kinds, "OmegaG$(g)")
        codon_param_vec = add_to_each_element(codon_param_vec, omegagrid)
    end
    if is_background
        push!(param_kinds, "OmegaBackground")
        codon_param_vec = add_to_each_element(codon_param_vec, background_omega_grid)
    end
    codon_param_vec

    num_sites = tree.message[1].sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l, num_sites)

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")


    # modify the tree 

    # compute recursion G1
    # compute recursion G2

    # How does the omega influence the propagation?

    # getnodelist and get a specific node from nodelist
    #function get_node_by_name(query_node::String, tree::FelNode)
    #    # This is a useful function that maybe should be in AbstractTreeNode.jl
    #    for node in getnodelist(tree)
    #        if node.name == query_node
    #            return node
    #        end
    #    end
    #end

    function check_purity_from_node_and_forward_in_time(tree)
        # collect all nodes under and see if they are the same group
        # if return true that means everything under the tree is pure (the tree node can be differnt though)
        node_groups = []
        for node in getnodelist(tree)
            if tree.name != node.name # we only check nodes under the input tree or subtree
                push!(node_groups, model_ind(node.name, tags))
            end
        end
        if length(unique(node_groups)) == 1
            return true
        else
            return false
        end
    end

    pure_clades = []
    function traverse_tree_to_check_for_pure_clades(pure_clades, tree)
        for child in tree.children
            if check_purity_from_node_and_forward_in_time(child)
                #print(child.name)
                #println(" PURE")
                # add this node, as everything under this node is pure
                push!(pure_clades, child)
            else
                #print(child.name)
                #println(" not pure") # keep recursion to se if we go forward in time if we can find pure clades
                traverse_tree_to_check_for_pure_clades(pure_clades, child)
            end
        end
        return pure_clades
    end

    time_elapsed = @elapsed subtree_tops = traverse_tree_to_check_for_pure_clades(pure_clades, tree)
    println("Patrick prune " * string(time_elapsed))

    cached_messages = Dict()
    for x in subtree_tops
        cached_messages[x] = Dict()
        parent = x.parent
        x.parent = nothing
        for cp in codon_param_vec
            alpha = cp[1]
            omegas = cp[2:end]

            relevant_omega = omegas[model_ind(x.children[1].name, tags)] # children[1].name is the same as children[2].name

            if haskey(cached_messages[x], (alpha, relevant_omega))
                continue
            end

            models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)
            felsenstein!(x, models)
            cached_messages[x][(alpha, relevant_omega)] = deepcopy(x.message)
        end
        x.parent = parent
    end

    for (row_ind, cp) in enumerate(codon_param_vec)
        alpha = cp[1]
        omegas = cp[2:end]
        tagged_models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)

        # added for subtree prune
        pure_children = []
        for x in subtree_tops
            push!(pure_children, x.children)
            x.message = cached_messages[x][(alpha, omegas[model_ind(x.children[1].name, tags)])]
            x.children = FelNode[]
        end

        felsenstein!(tree, tagged_models)
        #This combine!() is needed because the current site_LLs function applies to a partition
        #And after a felsenstein pass, you don't have the eq freqs factored in.
        #We could make a version of log_likelihood() that returns the partitions instead of just the sum
        combine!.(tree.message, tree.parent_message)

        log_con_lik_matrix[row_ind, :] .= MolecularEvolution.site_LLs(tree.message[1]) #Check that these grab the scaling constants as well!
        verbosity > 0 && if mod(row_ind, 500) == 1
            print(round(100 * row_ind / length(codon_param_vec)), "% ")
            flush(stdout)
        end

        for (i, x) in enumerate(subtree_tops)
            x.children = pure_children[i]
        end
    end

    verbosity > 0 && println()

    con_lik_matrix = zeros(size(log_con_lik_matrix))
    site_scalers = maximum(log_con_lik_matrix, dims=1)
    for i in 1:num_sites
        con_lik_matrix[:, i] .= exp.(log_con_lik_matrix[:, i] .- site_scalers[i])
    end

    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end

function difFUBAR_grid_final(tree, tags, GTRmat, F3x4_freqs, code; verbosity=1, foreground_grid=6, background_grid=4)

    cached_model = MG94_cacher(code)

    #This is the function that assigns models to branches
    #Sometimes there will be the same number of tags as omegas, but sometimes there will be one more omega.
    #Depending on whether there are any background branches (UNTESTED WITH NO BACKGROUND BRANCHES)
    function N_Omegas_model_func(tags, omega_vec, alpha, nuc_mat, F3x4, code)
        models = [cached_model(alpha, alpha * o, nuc_mat, F3x4, genetic_code=code) for o in omega_vec]
        return n::FelNode -> [models[model_ind(n.name, tags)]]
    end

    #Defines the grid used for inference.
    function gridsetup(lb, ub, num_below_one, trin, tr)
        step = (trin(1.0) - trin(lb)) / num_below_one
        return tr.(trin(lb):step:trin(ub))
    end
    tr(x) = 10^x - 0.05
    trinv(x) = log10(x + 0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    #is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree)]) > num_groups
    is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree) if !isroot(n)]) > num_groups

    tensor_dims = 1 + num_groups + is_background

    function add_to_each_element(vec_of_vec, elems)
        return [vcat(v, [e]) for v in vec_of_vec for e in elems]
    end

    codon_param_vec = [[a] for a in alphagrid]
    param_kinds = ["Alpha"]
    for g in 1:num_groups
        push!(param_kinds, "OmegaG$(g)")
        codon_param_vec = add_to_each_element(codon_param_vec, omegagrid)
    end
    if is_background
        push!(param_kinds, "OmegaBackground")
        codon_param_vec = add_to_each_element(codon_param_vec, background_omega_grid)
    end
    codon_param_vec

    num_sites = tree.message[1].sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l, num_sites)

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")

    function check_purity_from_node_and_forward_in_time(tree)
        # collect all nodes under and see if they are the same group
        # if return true that means everything under the tree is pure (the tree node can be differnt though)
        node_groups = []
        for node in getnodelist(tree)
            if tree.name != node.name # we only check nodes under the input tree or subtree
                push!(node_groups, model_ind(node.name, tags))
            end
        end
        if length(unique(node_groups)) == 1
            return true
        else
            return false
        end
    end

    pure_clades = []
    function traverse_tree_to_check_for_pure_clades(pure_clades, tree)
        for child in tree.children
            if check_purity_from_node_and_forward_in_time(child)
                # add this node, as everything under this node is pure
                push!(pure_clades, child)
            else
                traverse_tree_to_check_for_pure_clades(pure_clades, child)
            end
        end
        return pure_clades
    end

    time_elapsed = @elapsed subtree_tops = traverse_tree_to_check_for_pure_clades(pure_clades, tree)
    println("Patrick prune " * string(time_elapsed))

    cached_messages = Dict()
    cached_tag_inds = Dict()
    for x in subtree_tops
        cached_messages[x] = Dict()
        cached_tag_inds[x] = model_ind(x.children[1].name, tags)
        parent = x.parent
        x.parent = nothing
        for cp in codon_param_vec
            alpha = cp[1]
            omegas = cp[2:end]

            relevant_omega = omegas[cached_tag_inds[x]]

            if haskey(cached_messages[x], (alpha, relevant_omega))
                continue
            end

            models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)
            felsenstein!(x, models)
            cached_messages[x][(alpha, relevant_omega)] = deepcopy(x.message)
        end
        x.parent = parent
        x.children = FelNode[]
    end

    for (row_ind, cp) in enumerate(codon_param_vec)
        alpha = cp[1]
        omegas = cp[2:end]
        tagged_models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)

        # added for subtree prune 
        for x in subtree_tops
            x.message = cached_messages[x][(alpha, omegas[cached_tag_inds[x]])]
        end

        felsenstein!(tree, tagged_models)
        #This combine!() is needed because the current site_LLs function applies to a partition
        #And after a felsenstein pass, you don't have the eq freqs factored in.
        #We could make a version of log_likelihood() that returns the partitions instead of just the sum
        combine!.(tree.message, tree.parent_message)

        log_con_lik_matrix[row_ind, :] .= MolecularEvolution.site_LLs(tree.message[1]) #Check that these grab the scaling constants as well!
        verbosity > 0 && if mod(row_ind, 500) == 1
            print(round(100 * row_ind / length(codon_param_vec)), "% ")
            flush(stdout)
        end
    end

    verbosity > 0 && println()

    con_lik_matrix = zeros(size(log_con_lik_matrix))
    site_scalers = maximum(log_con_lik_matrix, dims=1)
    for i in 1:num_sites
        con_lik_matrix[:, i] .= exp.(log_con_lik_matrix[:, i] .- site_scalers[i])
    end

    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end


function difFUBAR_grid_pruned_3_return_subtree_tops(tree, tags, GTRmat, F3x4_freqs, code; verbosity=1, foreground_grid=6, background_grid=4)

    cached_model = MG94_cacher(code)

    #This is the function that assigns models to branches
    #Sometimes there will be the same number of tags as omegas, but sometimes there will be one more omega.
    #Depending on whether there are any background branches (UNTESTED WITH NO BACKGROUND BRANCHES)
    function N_Omegas_model_func(tags, omega_vec, alpha, nuc_mat, F3x4, code)
        models = [cached_model(alpha, alpha * o, nuc_mat, F3x4, genetic_code=code) for o in omega_vec]
        return n::FelNode -> [models[model_ind(n.name, tags)]]
    end

    #Defines the grid used for inference.
    function gridsetup(lb, ub, num_below_one, trin, tr)
        step = (trin(1.0) - trin(lb)) / num_below_one
        return tr.(trin(lb):step:trin(ub))
    end
    tr(x) = 10^x - 0.05
    trinv(x) = log10(x + 0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree)]) > num_groups
    tensor_dims = 1 + num_groups + is_background

    function add_to_each_element(vec_of_vec, elems)
        return [vcat(v, [e]) for v in vec_of_vec for e in elems]
    end

    codon_param_vec = [[a] for a in alphagrid]
    param_kinds = ["Alpha"]
    for g in 1:num_groups
        push!(param_kinds, "OmegaG$(g)")
        codon_param_vec = add_to_each_element(codon_param_vec, omegagrid)
    end
    if is_background
        push!(param_kinds, "OmegaBackground")
        codon_param_vec = add_to_each_element(codon_param_vec, background_omega_grid)
    end
    codon_param_vec

    num_sites = tree.message[1].sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l, num_sites)

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")


    # modify the tree 

    # compute recursion G1
    # compute recursion G2

    # How does the omega influence the propagation?

    # getnodelist and get a specific node from nodelist
    #function get_node_by_name(query_node::String, tree::FelNode)
    #    # This is a useful function that maybe should be in AbstractTreeNode.jl
    #    for node in getnodelist(tree)
    #        if node.name == query_node
    #            return node
    #        end
    #    end
    #end

    function check_purity_from_node_and_forward_in_time(tree)
        # collect all nodes under and see if they are the same group
        # if return true that means everything under the tree is pure (the tree node can be differnt though)
        node_groups = []
        for node in getnodelist(tree)
            if tree.name != node.name # we only check nodes under the input tree or subtree
                push!(node_groups, model_ind(node.name, tags))
            end
        end
        if length(unique(node_groups)) == 1
            return true
        else
            return false
        end
    end

    pure_clades = []
    function traverse_tree_to_check_for_pure_clades(pure_clades, tree)
        for child in tree.children
            if check_purity_from_node_and_forward_in_time(child)
                #print(child.name)
                #println(" PURE")
                # add this node, as everything under this node is pure
                push!(pure_clades, child)
            else
                #print(child.name)
                #println(" not pure") # keep recursion to se if we go forward in time if we can find pure clades
                traverse_tree_to_check_for_pure_clades(pure_clades, child)
            end
        end
        return pure_clades
    end

    time_elapsed = @elapsed subtree_tops = traverse_tree_to_check_for_pure_clades(pure_clades, tree)
    return subtree_tops
end

#foreground_grid and background_grid control the number of categories below 1.0
function difFUBAR_grid_pruned_1_return_subtree_tops(tree, tags, GTRmat, F3x4_freqs, code; verbosity=1, foreground_grid=6, background_grid=4)

    cached_model = MG94_cacher(code)

    #This is the function that assigns models to branches
    #Sometimes there will be the same number of tags as omegas, but sometimes there will be one more omega.
    #Depending on whether there are any background branches (UNTESTED WITH NO BACKGROUND BRANCHES)
    function N_Omegas_model_func(tags, omega_vec, alpha, nuc_mat, F3x4, code)
        models = [cached_model(alpha, alpha * o, nuc_mat, F3x4, genetic_code=code) for o in omega_vec]
        return n::FelNode -> [models[model_ind(n.name, tags)]]
    end

    #Defines the grid used for inference.
    function gridsetup(lb, ub, num_below_one, trin, tr)
        step = (trin(1.0) - trin(lb)) / num_below_one
        return tr.(trin(lb):step:trin(ub))
    end
    tr(x) = 10^x - 0.05
    trinv(x) = log10(x + 0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree)]) > num_groups
    tensor_dims = 1 + num_groups + is_background

    function add_to_each_element(vec_of_vec, elems)
        return [vcat(v, [e]) for v in vec_of_vec for e in elems]
    end

    codon_param_vec = [[a] for a in alphagrid]
    param_kinds = ["Alpha"]
    for g in 1:num_groups
        push!(param_kinds, "OmegaG$(g)")
        codon_param_vec = add_to_each_element(codon_param_vec, omegagrid)
    end
    if is_background
        push!(param_kinds, "OmegaBackground")
        codon_param_vec = add_to_each_element(codon_param_vec, background_omega_grid)
    end
    codon_param_vec

    num_sites = tree.message[1].sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l, num_sites)

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")


    # modify the tree 

    # compute recursion G1
    # compute recursion G2

    # How does the omega influence the propagation?

    pure_subclades = FelNode[]
    # Traverses the tree recursively with a dfs whilst pushing roots of pure subclades to list^
    function find_pure_subclades(node)
        # Get the index of the node's tag
        tag_ind_of_node = model_ind(node.name, tags)

        # If the node is a leaf, it's pure
        if length(node.children) == 0
            return true, tag_ind_of_node
        end

        children_is_pure_and_tag = []
        pure_children = FelNode[]

        for child in node.children
            child_is_pure, tag_ind = find_pure_subclades(child)
            if child_is_pure
                push!(pure_children, child)
            end
            push!(children_is_pure_and_tag, (child_is_pure, tag_ind))
        end

        # Get the index of the node's first child's tag
        tag_ind_of_first_child = children_is_pure_and_tag[1][2]

        # This is the case where the subclade starting at node is pure
        if all(x == (true, tag_ind_of_first_child) for x in children_is_pure_and_tag)
            if tag_ind_of_node != tag_ind_of_first_child
                # The purity is broken at this node
                push!(pure_subclades, node)
                return false, tag_ind_of_node
            end
            # The purity is not broken at this node
            return true, tag_ind_of_node
        end

        # This is the case where some child has mixed tags or the children are pure with regards to different tags
        for pure_child in pure_children
            if length(pure_child.children) == 0
                # We don't want to push leaves into pure_subclades
                continue
            end
            push!(pure_subclades, pure_child)
        end
        return false, tag_ind_of_node
    end

    time_elapsed = @elapsed find_pure_subclades(tree)
    return pure_subclades
end


function difFUBAR_sample(con_lik_matrix, iters; verbosity=1)
    verbosity > 0 && println("Step 4: Running Gibbs sampler to infer site categories.")
    alloc_grid, theta = LDA_gibbs_track_allocation_vec(con_lik_matrix, 0.1, iters=iters)
    return alloc_grid, theta
end

function difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=1, sites_to_plot=nothing, exports=true)
    grid_size, num_sites = size(alloc_grid)

    r(s) = round(s, digits=4)

    detected_sites = Int64[]
    group1_volumes = Vector{Float64}[]
    group2_volumes = Vector{Float64}[]
    alpha_volumes = Vector{Float64}[]
    detections = Vector{Float64}[] #legacy name - now includes all 4 "relevant" site posteriors
    param_means = Vector{Float64}[]

    ω1 = [c[2] for c in codon_param_vec]
    ω2 = [c[3] for c in codon_param_vec]
    alphas = [c[1] for c in codon_param_vec]
    ω1_greater_filt = ω1 .> ω2
    ω2_greater_filt = ω2 .> ω1
    ω1_pos_filt = ω1 .> 1.0
    ω2_pos_filt = ω2 .> 1.0

    verbosity > 0 && println("Step 5: Tabulating and plotting. Detected sites:")
    for site in 1:num_sites
        ω1_greater_posterior = sum(alloc_grid[ω1_greater_filt, site]) / sum(alloc_grid[:, site])
        ω2_greater_posterior = sum(alloc_grid[ω2_greater_filt, site]) / sum(alloc_grid[:, site])
        ω1_pos_posterior = sum(alloc_grid[ω1_pos_filt, site]) / sum(alloc_grid[:, site])
        ω2_pos_posterior = sum(alloc_grid[ω2_pos_filt, site]) / sum(alloc_grid[:, site])
        detecs = [ω1_greater_posterior, ω2_greater_posterior, ω1_pos_posterior, ω2_pos_posterior]

        site_counts_ω1 = collapse_counts(ω1, alloc_grid[:, site], cases=omegagrid)
        site_counts_ω2 = collapse_counts(ω2, alloc_grid[:, site], cases=omegagrid)
        site_counts_alphas = collapse_counts(alphas, alloc_grid[:, site], cases=alphagrid)

        mean_alpha = sum(site_counts_alphas .* alphagrid)
        mean_ω1 = sum(site_counts_ω1 .* omegagrid)
        mean_ω2 = sum(site_counts_ω2 .* omegagrid)

        push!(detections, detecs)
        push!(param_means, [mean_alpha, mean_ω1, mean_ω2])
        push!(group1_volumes, site_counts_ω1)
        push!(group2_volumes, site_counts_ω2)
        push!(alpha_volumes, site_counts_alphas)

        if maximum(detecs) > pos_thresh
            verbosity > 0 && print("Site $(site) - ")
            verbosity > 0 && print("P(ω1 > ω2):", ω1_greater_posterior)
            verbosity > 0 && print("; P(ω2 > ω1):", ω2_greater_posterior)
            verbosity > 0 && print("; P(ω1 > 1):", ω1_pos_posterior)
            verbosity > 0 && println("; P(ω2 > 1):", ω2_pos_posterior)
            push!(detected_sites, site)
        end
    end

    #Exporting site data
    df = DataFrame()
    df[!, "Codon Sites"] = [1:num_sites;]
    df[!, "P(ω1 > ω2)"] = [d[1] for d in detections]
    df[!, "P(ω2 > ω1)"] = [d[2] for d in detections]
    df[!, "P(ω1 > 1)"] = [d[3] for d in detections]
    df[!, "P(ω2 > 1)"] = [d[4] for d in detections]
    df[!, "mean(α)"] = [d[1] for d in param_means]
    df[!, "mean(ω1)"] = [d[2] for d in param_means]
    df[!, "mean(ω2)"] = [d[3] for d in param_means]

    verbosity > 0 && println("\nIf exports = true, writing results for all sites to CSV: " * analysis_name * "_posteriors.csv")
    exports && CSV.write(analysis_name * "_posteriors.csv", df)

    sites = [1:num_sites;]

    #Select the sites that will get plotted, in case you want to customize this.
    if isnothing(sites_to_plot)
        sites_to_plot = detected_sites
    end

    if length(sites_to_plot) == 0
        verbosity > 0 && println("No sites detected above threshold.")
    elseif exports
        verbosity > 0 && println("Plotting alpha and omega distributions. If exports = true, saved as " * analysis_name * "_violin_*.pdf")

        #Assumes alpha and omega grids are the same!? Currently enforced by args passed into difFUBAR_grid
        #Maybe this is ok
        grd = round.(omegagrid, digits=3)

        #Three plotting examples.
        #Plot the alphas for each flagged site

        lmargin = 7 + length(sites_to_plot) / 2
        ysize = 300 + 70 * length(sites[sites_to_plot])
        FUBAR_violin_plot(sites[sites_to_plot], alpha_volumes[sites_to_plot] .* 0.75, grd, tag="α", color="green", x_label="α")
        plot!(size=(400, ysize), grid=false, left_margin=(lmargin)mm, bottom_margin=10mm)

        savefig(analysis_name * "_violin_alpha.pdf")
        Plots.CURRENT_PLOT.nullableplot = nothing # PyPlots close()

        #Plot the G1 and G2 omegas
        FUBAR_violin_plot(sites[sites_to_plot], group1_volumes[sites_to_plot], grd, tag="ω1", color=tag_colors[1])
        FUBAR_violin_plot(sites[sites_to_plot], group2_volumes[sites_to_plot], grd, tag="ω2", color=tag_colors[2], x_label="ω")
        plot!(size=(400, ysize), grid=false, left_margin=(lmargin)mm, bottom_margin=10mm)

        savefig(analysis_name * "_violin_omegas.pdf")
        Plots.CURRENT_PLOT.nullableplot = nothing

        #Plot all three parameters, using the v_offset to separate the alphas from the omegas
        FUBAR_violin_plot(sites[sites_to_plot], group1_volumes[sites_to_plot] .* 0.5, grd, tag="ω1", color=tag_colors[1], v_offset=-0.1)
        FUBAR_violin_plot(sites[sites_to_plot], group2_volumes[sites_to_plot] .* 0.5, grd, tag="ω2", color=tag_colors[2], v_offset=-0.1)
        FUBAR_violin_plot(sites[sites_to_plot], alpha_volumes[sites_to_plot] .* 0.5, grd, tag="α", color="green", v_offset=0.1)
        plot!(size=(400, ysize), grid=false, left_margin=(lmargin)mm, bottom_margin=10mm)

        savefig(analysis_name * "_violin_all_params.pdf")
        Plots.CURRENT_PLOT.nullableplot = nothing

        #Coerce the violin plot function to also viz the "detection" posteriors.
        floored_detec = [clamp.((d .- 0.95) .* 20, 0.0, 1.0) for d in detections[sites_to_plot]]
        println(sites_to_plot)
        FUBAR_violin_plot(sites[sites_to_plot], [[f[1], 0.0, 0.0, 0.0] for f in floored_detec] .* 0.5,
            ["P(ω1>ω2)", "P(ω2>ω1)", "P(ω1>1)", "P(ω2>1)"], tag="P(ω1>ω2)", color=tag_colors[1],
            vertical_ind=nothing, plot_legend=false)
        FUBAR_violin_plot(sites[sites_to_plot], [[0.0, f[2], 0.0, 0.0] for f in floored_detec] .* 0.5,
            ["P(ω1>ω2)", "P(ω2>ω1)", "P(ω1>1)", "P(ω2>1)"], tag="P(ω2>ω1)", color=tag_colors[2],
            vertical_ind=nothing, plot_legend=false)
        FUBAR_violin_plot(sites[sites_to_plot], [[0.0, 0.0, f[3], 0.0] for f in floored_detec] .* 0.5,
            ["P(ω1>ω2)", "P(ω2>ω1)", "P(ω1>1)", "P(ω2>1)"], tag="P(ω1>1)", color=tag_colors[1],
            vertical_ind=nothing, plot_legend=false)
        FUBAR_violin_plot(sites[sites_to_plot], [[0.0, 0.0, 0.0, f[4]] for f in floored_detec] .* 0.5,
            ["P(ω1>ω2)", "P(ω2>ω1)", "P(ω1>1)", "P(ω2>1)"], tag="P(ω2>1)", color=tag_colors[2],
            vertical_ind=nothing, legend_ncol=2, x_label="", plot_legend=false)

        lmargin_detect = 12 + length(sites_to_plot) / 2

        plot!(size=(800, ysize), margins=1Plots.cm, legend=false, grid=false,
            ytickfont=18, bottom_margin=30mm, left_margin=(lmargin_detect)mm,
            xtickfont=18)
        println(length(sites_to_plot))

        savefig(analysis_name * "_detections.pdf")
        Plots.CURRENT_PLOT.nullableplot = nothing

    end

    if exports
        Plots.CURRENT_PLOT.nullableplot = nothing
        FUBAR_omega_plot(param_means, tag_colors, pos_thresh, detections, num_sites)
        plot!(size=(1.25 * length(sites), 300), margins=1Plots.cm, grid=false, legendfontsize=8)
        savefig(analysis_name * "_site_omega_means.pdf")

    end

    return df
end
export difFUBAR_tabulate

#Must return enough to re-calculate detections etc
function difFUBAR(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code)
    analysis_name = outpath
    tree, tags, tag_colors, analysis_name = difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit(seqnames, seqs, tree, generate_tag_stripper(tags), code, verbosity=verbosity)
    con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code,
        verbosity=verbosity, foreground_grid=6, background_grid=4)
    alloc_grid, theta = difFUBAR_sample(con_lik_matrix, iters, verbosity=verbosity)
    df = difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=verbosity, exports=exports)

    #Return df, (tuple of partial calculations needed to re-run tablulate)
    return df, (alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors)
end
export difFUBAR

#Must return enough to re-calculate detections etc
function difFUBAR_prune_max(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code)
    analysis_name = outpath
    tree, tags, tag_colors, analysis_name = difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit(seqnames, seqs, tree, generate_tag_stripper(tags), code, verbosity=verbosity)
    con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = difFUBAR_grid_pruned_1(tree, tags, GTRmat, F3x4_freqs, code,
        verbosity=verbosity, foreground_grid=6, background_grid=4)
    alloc_grid, theta = difFUBAR_sample(con_lik_matrix, iters, verbosity=verbosity)
    df = difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=verbosity, exports=exports)

    #Return df, (tuple of partial calculations needed to re-run tablulate)
    return df, (alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors)
end
export difFUBAR_prune_max


#Must return enough to re-calculate detections etc
function difFUBAR_prune_patrick(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code)
    analysis_name = outpath
    tree, tags, tag_colors, analysis_name = difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit(seqnames, seqs, tree, generate_tag_stripper(tags), code, verbosity=verbosity)
    con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = difFUBAR_grid_pruned_2(tree, tags, GTRmat, F3x4_freqs, code,
        verbosity=verbosity, foreground_grid=6, background_grid=4)
    alloc_grid, theta = difFUBAR_sample(con_lik_matrix, iters, verbosity=verbosity)
    df = difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=verbosity, exports=exports)

    #Return df, (tuple of partial calculations needed to re-run tablulate)
    return df, (alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors)
end
export difFUBAR_prune_patrick

function difFUBAR_prune_patrick_max(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code)
    analysis_name = outpath
    tree, tags, tag_colors, analysis_name = difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit(seqnames, seqs, tree, generate_tag_stripper(tags), code, verbosity=verbosity)
    con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = difFUBAR_grid_pruned_3(tree, tags, GTRmat, F3x4_freqs, code,
        verbosity=verbosity, foreground_grid=6, background_grid=4)
    alloc_grid, theta = difFUBAR_sample(con_lik_matrix, iters, verbosity=verbosity)
    df = difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=verbosity, exports=exports)

    #Return df, (tuple of partial calculations needed to re-run tablulate)
    return df, (alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors)
end
export difFUBAR_prune_patrick_max

function difFUBAR_prune_patrick_max_child(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code)
    analysis_name = outpath
    tree, tags, tag_colors, analysis_name = difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit(seqnames, seqs, tree, generate_tag_stripper(tags), code, verbosity=verbosity)
    con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = difFUBAR_grid_pruned_4(tree, tags, GTRmat, F3x4_freqs, code,
        verbosity=verbosity, foreground_grid=6, background_grid=4)
    alloc_grid, theta = difFUBAR_sample(con_lik_matrix, iters, verbosity=verbosity)
    df = difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=verbosity, exports=exports)

    #Return df, (tuple of partial calculations needed to re-run tablulate)
    return df, (alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors)
end
export difFUBAR_prune_patrick_max_child

function difFUBAR_prune_final(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code)
    analysis_name = outpath
    tree, tags, tag_colors, analysis_name = difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit(seqnames, seqs, tree, generate_tag_stripper(tags), code, verbosity=verbosity)
    con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = difFUBAR_grid_final(tree, tags, GTRmat, F3x4_freqs, code,
        verbosity=verbosity, foreground_grid=6, background_grid=4)
    alloc_grid, theta = difFUBAR_sample(con_lik_matrix, iters, verbosity=verbosity)
    df = difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=verbosity, exports=exports)

    #Return df, (tuple of partial calculations needed to re-run tablulate)
    return df, (alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors)
end
export difFUBAR_prune_final

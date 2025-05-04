"""
This file contains different versions of the difFUBAR_grid algorithm.
Including: the default difFUBAR_grid function as well as different speedup tricks related to it, i.e.,
combinations of implementations using parallelization and memoization. Also, it contains helper functions for these.
"""

"""
# Constructor
    difFUBARBaseline()

# Description
Use the trivial implementation of the grid likelihood computations, i.e. 1 thread without sub-tree likelihood caching.

See also: `difFUBARParallel`, `difFUBARTreesurgery`, `difFUBARTreesurgeryAndParallel`.
"""
struct difFUBARBaseline <: difFUBARGrid end
"""
# Constructor
    difFUBARParallel()

# Description
Extend the baseline version by parallelizing the grid calculations. Requires julia to be launched with the `t` switch.
Using `t` computational threads, where `t` is sufficiently small, memory complexity is usually O(t) and time complexity O(1/t).
Empirical tests suggests that `t` should not be higher than the machine's total CPU threads and usually not higher than half of it's total threads.

See also: `difFUBARBaseline`, `difFUBARTreesurgery`, `difFUBARTreesurgeryAndParallel`.
"""
struct difFUBARParallel <: difFUBARGrid end
"""
# Constructor
    difFUBARTreesurgery()

# Description
Use sub-tree likelihood caching described in the "Methods" section of the difFUBAR paper.
Use more memory than the baseline version but be significantly faster,
if purity is high.

See also: `difFUBARBaseline`, `difFUBARParallel`, `difFUBARTreesurgeryAndParallel`.
"""
struct difFUBARTreesurgery <: difFUBARGrid end
"""
# Constructor
    difFUBARTreesurgeryAndParallel()

# Description
Use parallelization and sub-tree likelihood caching. The most performant version in most cases. Use more memory than other versions.

See also: `difFUBARBaseline`, `difFUBARTreesurgery`, `difFUBARParallel`.
"""
struct difFUBARTreesurgeryAndParallel <: difFUBARGrid end

function add_to_each_element(vec_of_vec, elems)
    return [vcat(v,[e]) for v in vec_of_vec for e in elems]
end

function generate_alpha_and_single_omega_grids(alphagrid, omegagrid, background_omega_grid, is_background)
    alpha_and_single_omega_grids = Dict()
    alphagrid_vectorized = [[a] for a in alphagrid]
    alpha_and_single_omega_grids["Omega"] = add_to_each_element(alphagrid_vectorized,omegagrid)
    if is_background
        alpha_and_single_omega_grids["OmegaBackground"] = add_to_each_element(alphagrid_vectorized,background_omega_grid)
    end
    return alpha_and_single_omega_grids
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

#Precalculates the DiagonalizedCTMC for all unique alpha and beta pairs and puts them in the cached model. Used in the parallel versions.
function precalculate_models!(cached_model, alphagrid, omegagrid, background_omega_grid, is_background, GTRmat, F3x4_freqs)
    if is_background
        unique_omega_grid = unique(vcat(omegagrid, background_omega_grid))
    else
        unique_omega_grid = omegagrid
    end
    for alpha in alphagrid
        for omega in unique_omega_grid
            cached_model(alpha, alpha*omega, GTRmat, F3x4_freqs)
        end
    end
end
"""
    getpuresubclades(tree::FelNode, tags::Vector{String})

- Should usually be called on the root of the tree. Traverses the tree iteratively with a depth-first search to find roots of pure subclades, presuming that nodenames have been trailed with tags. Returns a Vector{FelNode} with root-nodes of the pure subclades.
"""
function getpuresubclades(tree::FelNode, tags::Vector{String})
    pure_subclades = Vector{FelNode}()
    node_stack = [(tree, model_ind(tree.name, tags), true)]
    result_stack = Vector{Tuple{Bool, Int64}}() #(is_pure, tag_ind)
    while !isempty(node_stack)
        node, tag_ind_of_node, first_time = pop!(node_stack)
        if !isleafnode(node)
            if first_time
                push!(node_stack, (node, tag_ind_of_node, false))
                for child in node.children
                    push!(node_stack, (child, model_ind(child.name, tags), true))
                end
            end
            if !first_time
                children_are_pure = Vector{Bool}()
                children_tag_inds = Vector{Int64}()

                for _ = 1:length(node.children)
                    #children_are_pure etc. maps directly onto node.children since we have two pop! operations (reverses order) that cancel out
                    child_is_pure, tag_ind = pop!(result_stack)
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
                        push!(result_stack, (false, tag_ind_of_node))
                        continue
                    end
                    # The purity is not broken at this node
                    push!(result_stack, (true, tag_ind_of_node))
                    continue
                end

                # This is the case where some child has mixed tags or the children are pure with regards to different tags
                for (child_is_pure, child) in zip(children_are_pure, node.children)
                    if !child_is_pure || isleafnode(child)
                        # We don't want to push leaves into pure_subclades
                        continue
                    end
                    push!(pure_subclades, child)
                end
                push!(result_stack, (false, tag_ind_of_node))
            end
        else
            # If the node is a leaf, it's pure
            push!(result_stack, (true, tag_ind_of_node))
        end
    end
    return pure_subclades
end

#Calculates the ratio of nodes that are in a pure clade to total nodes in the tree (1st return value)
#Calculates the number of pure omega and background omega clades (2nd and 3rd return values)
function get_purity_info(tree, tags, num_groups)
    pure_subclades = getpuresubclades(tree, tags)
    c = 0
    for x in pure_subclades
        #The root of the pure clade is not counted
        c += length(getnodelist(x)) - 1
    end
    num_omega_clades = count(x -> model_ind(x.children[1].name, tags) <= num_groups, pure_subclades)
    num_background_omega_clades = length(pure_subclades) - num_omega_clades
    return c / length(getnodelist(tree)), num_omega_clades, num_background_omega_clades
end

#Estimates the extra memory that the tree-surgery version would use by caching messages
function extra_mem_from_caching(num_sites, alphagrid, omegagrid, background_omega_grid, num_omega_clades, num_background_omega_clades, code)
    num_alphas, num_omegas, num_background_omegas = map(length, (alphagrid, omegagrid, background_omega_grid))
    F64s_in_CodonPartition = num_sites * (length(code.sense_codons) + 1) # the '+ 1' is for .scaling
    total_bytes_from_caching = 8 * F64s_in_CodonPartition * num_alphas * (num_omegas * num_omega_clades + num_background_omegas * num_background_omega_clades)
    return total_bytes_from_caching
end

function choose_grid_and_nthreads(tree, tags, num_groups, num_sites, alphagrid, omegagrid, background_omega_grid, code)
    purity_ratio, num_omega_clades, num_background_omega_clades = get_purity_info(tree, tags, num_groups)
    extra_mem = extra_mem_from_caching(num_sites, alphagrid, omegagrid, background_omega_grid, num_omega_clades, num_background_omega_clades, code)
    
    tree_size = Base.summarysize(tree) #Estimates (often overestimates) memory usage of one tree
    if Sys.isapple()
        free_memory = Sys.total_memory() - (parse(Int, chomp(read(`ps -o rss= -p $(getpid())`, String))) * 1024) #Assumes that other processes aren't using a lot of RAM
        optimal_nthreads_for_system = Sys.CPU_THREADS #Only counts performance cores
    else
        free_memory = Sys.free_memory()
        optimal_nthreads_for_system = Sys.CPU_THREADS ÷ 2
    end
    max_nthreads = Int(min(free_memory ÷ tree_size + 1, Threads.nthreads(), optimal_nthreads_for_system))
    tree_surgery_is_considered = extra_mem < free_memory && purity_ratio > 0.1 #Trying not to introduce too much overhead if the speedup isn't significant
    parallelization_is_considered = max_nthreads > 1

    if !(tree_surgery_is_considered || parallelization_is_considered)
        return difFUBARBaseline(), 1
    elseif !tree_surgery_is_considered
        return difFUBARParallel(), max_nthreads
    elseif !parallelization_is_considered
        return difFUBARTreesurgery(), 1
    end
    #From now on, we know that both tree-surgery and parallelization are considered
    if extra_mem + tree_size * max_nthreads < free_memory
        return difFUBARTreesurgeryAndParallel(), max_nthreads
    end
    if purity_ratio < 0.5 #In this case, maximum speedup from tree-surgery would be 2x, but we know that max_nthreads > 1
        return difFUBARParallel(), max_nthreads
    end
    max_nthreads_for_treesurgery_and_parallel = Int((free_memory - extra_mem) ÷ tree_size) + 1
    if max_nthreads_for_treesurgery_and_parallel > 1
        return difFUBARTreesurgeryAndParallel(), max_nthreads_for_treesurgery_and_parallel
    end
    #Now we must choose between tree surgery and parallelization
    if 1 - purity_ratio < 1 / max_nthreads #Estimates and compares speedup coefficients
        return difFUBAR_grid_treesurgery, 1
    else
        return difFUBARParallel(), max_nthreads
    end
end

#Defines the grid used for inference.
function gridsetup(lb, ub, num_below_one, trin, tr)
    step = (trin(1.0) - trin(lb))/num_below_one
    return tr.(trin(lb):step:trin(ub))
end

#Initializes variables common to all grid versions of difFUBAR
function gridprep(tree, tags; verbosity = 1, foreground_grid = 6, background_grid = 4)
    tr(x) = 10^x-0.05
    trinv(x) =  log10(x+0.05)
    alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr); 
    omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
    background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
    length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

    num_groups = length(tags)
    is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree) if !MolecularEvolution.isroot(n)]) > num_groups
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
    
    num_sites = tree.parent_message[1].partition.sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l,num_sites);
    return log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites
end

#Runs felsenstein! on a subgrid (an enumerated codon_param_vec chunk) and puts the results in log_con_lik_matrix. 
#Used in the parallel version.
function do_subgrid!(tree::FelNode, cached_model, cpv_chunk::Vector{Tuple{Int64, Vector{Float64}}}, tags::Vector{String}, GTRmat, F3x4_freqs, code, log_con_lik_matrix)
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

#Same as above but adapted for the version that does both the tree-surgery memoization and parallelization.
function do_subgrid!(tree::FelNode, cached_model, cpv_chunk::Vector{Tuple{Int64, Vector{Float64}}}, idx::Int64, pure_subclades::Vector{FelNode}, nodelists::Vector{Vector{FelNode}}, cached_messages, cached_tag_inds, tags::Vector{String}, GTRmat, F3x4_freqs, code, log_con_lik_matrix)
    # Note that cpv_chunk is already enumerated
    for (row_ind,cp) in cpv_chunk
        alpha = cp[1]
        omegas = cp[2:end]
        tagged_models = N_Omegas_model_func(cached_model, tags,omegas,alpha,GTRmat,F3x4_freqs, code)

        for x in pure_subclades
            nodeindex = x.nodeindex
            # Get the local equivalent node to x
            y = nodelists[idx][nodeindex]
            y.message[1].partition = cached_messages[nodeindex][(alpha, omegas[cached_tag_inds[nodeindex]])]
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

#Precalculates and returns messages for a pure subclade for different alpha and omega pairs. Aso is short for "alpha and single omega". 
#Used in the version that does both the tree-surgery memoization and parallelization.
function get_messages_for_aso_pairs(pure_subclade::FelNode, cached_model, aso_chunk::SubArray{Vector{Float64}}, GTRmat, F3x4_freqs, code)
    cached_messages_x = Dict()
    for (alpha, omega) in aso_chunk
        model = Omega_model_func(cached_model,omega,alpha,GTRmat,F3x4_freqs,code)
        felsenstein!(pure_subclade, model)
        cached_messages_x[(alpha, omega)] = copy_partition(pure_subclade.message[1].partition)
        MolecularEvolution.safe_release_partition!(pure_subclade.message[1])
    end
    return cached_messages_x
end

function difFUBAR_grid(version::difFUBARBaseline, tree, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites, nthreads; verbosity = 1, foreground_grid = 6, background_grid = 4)
    cached_model = MG94_cacher(code)
    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.parent_message[1].partition.sites) conditional likelihood values (the slowest step). Currently on:")

    for (row_ind,cp) in enumerate(codon_param_vec)
        alpha = cp[1]
        omegas = cp[2:end]
        tagged_models = N_Omegas_model_func(cached_model, tags,omegas,alpha,GTRmat,F3x4_freqs, code)

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
    alpha = codon_param_vec[1][1]
    o = codon_param_vec[1][2]

    verbosity > 0 && println()

    con_lik_matrix = zeros(size(log_con_lik_matrix));
    site_scalers = maximum(log_con_lik_matrix, dims = 1);
    for i in 1:num_sites
        con_lik_matrix[:,i] .= exp.(log_con_lik_matrix[:,i] .- site_scalers[i])
    end

    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end

function difFUBAR_grid(version::difFUBARParallel, tree, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites, nthreads; verbosity = 1, foreground_grid = 6, background_grid = 4)
    BLAS_num_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1) #Otherwise, BLAS threads inhibit Julia threads
    cached_model = MG94_cacher(code)
    precalculate_models!(cached_model, alphagrid, omegagrid, background_omega_grid, is_background, GTRmat, F3x4_freqs)
    trees = [tree, [deepcopy(tree) for _ = 1:(nthreads - 1)]...]
    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.parent_message[1].partition.sites) conditional likelihood values (the slowest step). Currently on:")

    cpv_chunks = Iterators.partition(enumerate(codon_param_vec), max(1, ceil(Int, length(codon_param_vec) / nthreads)))
    tasks = []
    for (i, cpv_chunk) in enumerate(cpv_chunks)
        # Spawn the task and add it to the array
        task = Threads.@spawn do_subgrid!(trees[i], cached_model, cpv_chunk, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix)
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
    BLAS.set_num_threads(BLAS_num_threads)
    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end

function difFUBAR_grid(version::difFUBARTreesurgery, tree, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites, nthreads; verbosity = 1, foreground_grid = 6, background_grid = 4)
    MolecularEvolution.set_node_indices!(tree)
    cached_model = MG94_cacher(code)

    pure_subclades = getpuresubclades(tree, tags)

    if length(pure_subclades) > 0
        alpha_and_single_omega_grids = generate_alpha_and_single_omega_grids(alphagrid, omegagrid, background_omega_grid, is_background)
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
        cached_messages[nodeindex] = Dict()
        parent = x.parent
        x.parent = nothing
        for (alpha, omega) in alpha_and_single_omega_grid
            model = Omega_model_func(cached_model,omega,alpha,GTRmat,F3x4_freqs,code)
            felsenstein!(x, model)
            cached_messages[nodeindex][(alpha, omega)] = copy_partition(x.message[1].partition)
            MolecularEvolution.safe_release_partition!(x.message[1])
        end
        x.message[1].static = true # Later on, we don't want the partition field to be recycled
        x.parent = parent
        x.children = FelNode[]
    end

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.parent_message[1].partition.sites) conditional likelihood values (the slowest step). Currently on:")

    for (row_ind,cp) in enumerate(codon_param_vec)
        alpha = cp[1]
        omegas = cp[2:end]
        tagged_models = N_Omegas_model_func(cached_model,tags,omegas,alpha,GTRmat,F3x4_freqs, code)

        for x in pure_subclades
            nodeindex = x.nodeindex
            x.message[1].partition = cached_messages[nodeindex][(alpha, omegas[cached_tag_inds[nodeindex]])]
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

function difFUBAR_grid(version::difFUBARTreesurgeryAndParallel, tree, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites, nthreads; verbosity = 1, foreground_grid = 6, background_grid = 4)
    BLAS_num_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1) #Otherwise, BLAS threads inhibit Julia threads
    MolecularEvolution.set_node_indices!(tree)
    trees = [tree, [deepcopy(tree) for _ = 1:(nthreads - 1)]...]
    cached_model = MG94_cacher(code)
    precalculate_models!(cached_model, alphagrid, omegagrid, background_omega_grid, is_background, GTRmat, F3x4_freqs)
    nodelists = [getnodelist(tree) for tree in trees]
    
    pure_subclades = getpuresubclades(tree, tags)

    if length(pure_subclades) > 0
        alpha_and_single_omega_grids = generate_alpha_and_single_omega_grids(alphagrid, omegagrid, background_omega_grid, is_background)
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

        aso_chunks = Iterators.partition(alpha_and_single_omega_grid, max(1, ceil(Int, length(alpha_and_single_omega_grid) / nthreads)))
        tasks = []
        for (i, aso_chunk) in enumerate(aso_chunks)
            # Spawn the task and add it to the array
            task = Threads.@spawn get_messages_for_aso_pairs(nodelists[i][nodeindex], cached_model, aso_chunk, GTRmat, F3x4_freqs, code)
            push!(tasks, task)
        end

        # Wait for all tasks to finish and collect their return values
        chunks_of_cached_messages = [fetch(task) for task in tasks];

        # Merge all the returned Dicts, and put it in the big Dict of messages
        cached_messages[nodeindex] = merge(chunks_of_cached_messages...)

        # Make sure that all of these partition fields are not recycled (this must be done after fetch)
        for j in eachindex(nodelists)
            # Index into this pure subclade in every tree
            nodelists[j][nodeindex].message[1].static = true
        end

        for (nodelist, parent) in zip(nodelists, parents)
            nodelist[nodeindex].parent = parent
            nodelist[nodeindex].children = FelNode[]
        end
    end

    GC.gc()

    num_sites = tree.parent_message[1].partition.sites
    l = length(codon_param_vec)
    log_con_lik_matrix = zeros(l,num_sites);

    verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.parent_message[1].partition.sites) conditional likelihood values (the slowest step). Currently on:")

    cpv_chunks = Iterators.partition(enumerate(codon_param_vec), max(1, ceil(Int, length(codon_param_vec) / nthreads)))
    tasks = []
    for (i, cpv_chunk) in enumerate(cpv_chunks)
        # Spawn the task and add it to the array
        task = Threads.@spawn do_subgrid!(trees[i], cached_model, cpv_chunk, i, pure_subclades, nodelists, cached_messages, cached_tag_inds, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix)
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
    BLAS.set_num_threads(BLAS_num_threads)
    return con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds
end
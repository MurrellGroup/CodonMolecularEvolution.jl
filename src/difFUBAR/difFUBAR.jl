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

const DIFFUBAR_TAG_COLORS = ["#ff0000", "#1900ff"] # [red, blue]

#Here, untagged comes after the tags
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

"""


"""
function difFUBAR_init(outpath_and_file_prefix, treestring, tags; tag_colors=DIFFUBAR_TAG_COLORS[sortperm(tags)], verbosity=1, exports=true, strip_tags_from_name=generate_tag_stripper(tags), disable_binarize=false, ladderize_tree = false)

    #Create the export directory, if required
    analysis_name = outpath_and_file_prefix
    splt = splitpath(analysis_name)[1:end-1]
    if length(splt) > 0
        exports && mkpath(joinpath(splt))
    end

    tree = gettreefromnewick(treestring, FelNode, disable_binarize=disable_binarize)

    #Need to think of consequences of eg. binarizing when there are tags.
    #We'll need the non-zero branch lengths to inherit the names/tags, for example.
    #REQUIRED TEST: CHECK NODE NAMES AFTER BINARIZATION
    #MolecularEvolution.binarize!(tree) #Check if this is required for trees with ternary branching?
    
    if ladderize_tree
        MolecularEvolution.ladderize!(tree)
    end

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
        plot_tagged_phylo_tree(PlotsExtDummy(), tree, tag_colors, tags, analysis_name)
    end

    #Tags and tag colors are now ordered, and tag_colors includes the untagged category
    return tree, tags, tag_colors, analysis_name
end

function difFUBAR_global_fit(seqnames, seqs, tree, leaf_name_transform, code; verbosity=1, optimize_branch_lengths=false)

    verbosity > 0 && println("Step 2: Optimizing global codon model parameters.")

    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = optimize_MG94_F3x4(seqnames, seqs, tree, leaf_name_transform=leaf_name_transform, genetic_code=code)

    ######
    #optionally polish branch lengths and topology
    ######

    return tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs
end

function difFUBAR_global_fit_2steps(seqnames, seqs, tree, leaf_name_transform, code; verbosity=1, optimize_branch_lengths=false)

    verbosity > 0 && println("Step 2: Optimizing global codon model parameters.")

    tree, nuc_mu, nuc_pi = optimize_nuc_mu(seqnames, seqs, tree, leaf_name_transform=leaf_name_transform, genetic_code=code, optimize_branch_lengths=optimize_branch_lengths)

    #Optionally polish branch lengths
    if optimize_branch_lengths == true
        tree_polish!(tree, GeneralCTMC(reversibleQ(nuc_mu, nuc_pi)), verbose=verbosity, topology=false)
        #Detect if all branchlengths are zero or all branchlengths are the same
    elseif optimize_branch_lengths == "detect"
        branchlengths = [x.branchlength for x in getnodelist(tree)]
        if all(x -> x == 0, branchlengths)
            @warn "All branchlengths are zero"
        elseif length(unique(branchlengths)) == 1
            @warn "All branchlengths are the same"
        end
    end

    GTRmat = reversibleQ(nuc_mu, ones(4))
    tree, alpha, beta, F3x4_freqs, eq_freqs = optimize_codon_alpha_and_beta(seqnames, seqs, tree, GTRmat, leaf_name_transform=leaf_name_transform, genetic_code=code)
    
    rescale_branchlengths!(tree, alpha) #rescale such that the ML value of alpha is 1

    ######
    #optionally polish branch lengths and topology
    ######

    return tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs
end

#foreground_grid and background_grid control the number of categories below 1.0
function difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code; verbosity=1, foreground_grid=6, background_grid=4, version::Union{difFUBARGrid,Nothing}=nothing, t=0)

    log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites = gridprep(tree, tags;
        verbosity=verbosity,
        foreground_grid=foreground_grid,
        background_grid=background_grid
    )
    heuristic_pick, nthreads = choose_grid_and_nthreads(tree, tags, num_groups, num_sites, alphagrid, omegagrid, background_omega_grid, code)
    t < 1 && (t = nthreads)
    isnothing(version) && (version = heuristic_pick)
    verbosity > 1 && println("\n$(typeof(version)), with $(t) parallel threads will be used for the grid likelihood calculations\n") #Higher level of verbosity
    return difFUBAR_grid(version, tree, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites, t;
        verbosity=verbosity,
        foreground_grid=foreground_grid,
        background_grid=background_grid
    )
end

function difFUBAR_sample(con_lik_matrix, iters; verbosity=1)
    verbosity > 0 && println("Step 4: Running Gibbs sampler to infer site categories.")
    alloc_grid, theta = LDA_gibbs_track_allocation_vec(con_lik_matrix, 0.1, iters=iters)
    return alloc_grid, theta
end

function difFUBAR_bayesian_postprocessing(pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid; tag_colors=DIFFUBAR_TAG_COLORS, verbosity=1)
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
    
    # Return everything needed by both tabulate and plot functions
    return detections, param_means, detected_sites, group1_volumes, group2_volumes, alpha_volumes, num_sites
end


function difFUBAR_tabulate(analysis_name, detections, param_means, num_sites; tag_colors=DIFFUBAR_TAG_COLORS, verbosity=1, exports=true)

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
    
    return df
end

export difFUBAR_tabulate_and_plot
function difFUBAR_tabulate_and_plot(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid; tag_colors=DIFFUBAR_TAG_COLORS, verbosity=1, exports=true)
    # Process the data and get all needed values
    detections, param_means, detected_sites, group1_volumes, group2_volumes, alpha_volumes, num_sites = 
        difFUBAR_bayesian_postprocessing(pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid; 
                                        tag_colors=tag_colors, verbosity=verbosity)

    # Pass appropriate values to each function
    df = difFUBAR_tabulate(analysis_name, detections, param_means, num_sites; 
                          tag_colors=tag_colors, verbosity=verbosity, exports=exports)
    
    # Make sure to pass all required values to plot_results
    difFUBAR_plot_results(PlotsExtDummy(), analysis_name, pos_thresh, detections, param_means, num_sites, omegagrid,
                         detected_sites, group1_volumes, group2_volumes, alpha_volumes;
                         tag_colors=tag_colors, verbosity=verbosity, exports=exports)
    
    return df
end

#Must return enough to re-calculate detections etc
export difFUBAR
"""
    difFUBAR(seqnames, seqs, treestring, tags, outpath; <keyword arguments>)

Takes a tagged phylogeny and an alignment as input and performs difFUBAR analysis.

# Arguments
- `seqnames`: vector of untagged sequence names.
- `seqs`: vector of aligned sequences, corresponding to `seqnames`.
- `treestring`: a tagged newick tree string.
- `tags`: vector of tag signatures.
- `outpath`: export directory.
- `tag_colors=DIFFUBAR_TAG_COLORS[sortperm(tags)]`: vector of tag colors (hex format). The default option is consistent with the difFUBAR paper (Foreground 1: red, Foreground 2: blue).
- `pos_thresh=0.95`: threshold of significance for the posteriors.
- `iters=2500`: iterations used in the Gibbs sampler.
- `verbosity=1`: as verbosity increases, prints are added accumulatively. 
    - 0 - no prints
    - 1 - show current step and where output files are exported
    - 2 - show the chosen `difFUBAR_grid` version and amount of parallel threads.
- `exports=true`: if true, output files are exported.
- `code=MolecularEvolution.universal_code`: genetic code used for the analysis.
- `optimize_branch_lengths=false`: if true, the branch lengths of the phylogenetic tree are optimized.
- `version::Union{difFUBARGrid, Nothing}=nothing`: explicitly choose the version of `difFUBAR_grid` to use. If `nothing`, the version is heuristically chosen based on the available RAM and Julia threads.
- `t=0`: explicitly choose the amount of Julia threads to use. If `0`, the degree of parallelization is heuristically chosen based on the available RAM and Julia threads.

!!! note
    Julia starts up with a single thread of execution, by default. See [Starting Julia with multiple threads](https://docs.julialang.org/en/v1/manual/multi-threading/#Starting-Julia-with-multiple-threads).
"""
function difFUBAR(seqnames, seqs, treestring, tags, outpath; tag_colors=DIFFUBAR_TAG_COLORS[sortperm(tags)], pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code, optimize_branch_lengths=false, version::Union{difFUBARGrid,Nothing}=nothing, t=0)
    analysis_name = outpath
    tree, tags, tag_colors, analysis_name = difFUBAR_init(analysis_name, treestring, tags, tag_colors=tag_colors, exports=exports, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit_2steps(seqnames, seqs, tree, generate_tag_stripper(tags), code, verbosity=verbosity, optimize_branch_lengths=optimize_branch_lengths)
    con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code,
        verbosity=verbosity, foreground_grid=6, background_grid=4, version=version, t=t)
    alloc_grid, theta = difFUBAR_sample(con_lik_matrix, iters, verbosity=verbosity)
    df = difFUBAR_tabulate_and_plot(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid; tag_colors=tag_colors, verbosity=verbosity, exports=exports)

    #Return df, (tuple of partial calculations needed to re-run tablulate)
    return df, (alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors)
end



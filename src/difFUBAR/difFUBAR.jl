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

    bar!([-10], [1], bottom=[1000], color=color, alpha=alpha, label=tag, linewidth=0, bar_edges=false, linealpha=0.0)
    bar!(yticks=(ypos, ["       " for _ in sites]), xticks=((1:length(omegagrid)), ["       " for _ in 1:length(omegagrid)]), xrotation=90)
    annotate!(length(omegagrid)/2, -length(sites)/2-(2.0+(length(sites)/500)), Plots.text(x_label, "black", :center, 10))
    annotate!(-4.5, -(length(sites)+1)/4, Plots.text(y_label, "black", :center, 10, rotation=90))

    bar!(ylim=(minimum(ypos) - 0.5, maximum(ypos) + 0.5), xlim = (0, length(omegagrid) + 1))
    if plot_legend
        plot!(
            legend=(0.5, 1+1.5/(50+length(sites))),
            legendcolumns=legend_ncol,
            shadow=true, fancybox=true,
        )
    end
    for i in 1:length(sites)
        annotate!(-0.5, ypos[i], Plots.text("$(sites[i])", "black", :right, 9))
    end
    for i in 1:length(omegagrid) 
        annotate!(i, -length(sites)/2-0.55 -length(sites)/3000, Plots.text("$(omegagrid[i])", "black", :right, 9, rotation=90))
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
        #TODO: update plots in docs
        phylo_tree = get_phylo_tree(tree)
        tagging = [tag_colors[model_ind(n, tags)] for n in nodenameiter(phylo_tree)]
        for node in nodeiter(phylo_tree)
            renamenode!(phylo_tree, node, strip_tags_from_name(node.name))
        end
        #Warnings regarding marker- and linecolor also appear in the Phylo.jl docs example
        #Note: sometimes long leafnames are truncated/not visible in the plot
        pl = plot(phylo_tree,
            showtips = true, tipfont = 6, markercolor = tagging, linecolor = tagging, markerstrokewidth = 0, size = (600, (120 + length(getleaflist(tree)) * 8)))
        savefig_tweakSVG(analysis_name * "_tagged_input_tree.svg", pl)
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

export difFUBAR_tabulate
"""
    difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid; tag_colors=DIFFUBAR_TAG_COLORS, verbosity=1, sites_to_plot=nothing, exports=true)

Takes the output of `difFUBAR` and tabulates the results. Returns a DataFrame of the detected sites.
This function enables you to use the results of `difFUBAR` to tabulate the results with a different threshold.

# Arguments
- `analysis_name`: where to export the results.
- `pos_thresh`: threshold of significance for the posteriors.
- `alloc_grid`: contains the result of the Gibbs sampler.
- `codon_param_vec`: vector of codon parameters from difFUBAR.
- `alphagrid`: grid of alpha values.
- `omegagrid`: grid of omega values.
- `tag_colors`: colors of the tags.
"""
function difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid; tag_colors=DIFFUBAR_TAG_COLORS, verbosity=1, sites_to_plot=nothing, exports=true)
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
        #FUBAR_violin_plot(sites[sites_to_plot], alpha_volumes[sites_to_plot] .* 0.75, grd, tag="α", color="green", x_label="α")
        Plots.CURRENT_PLOT.nullableplot = nothing # PyPlots close()
        FUBAR_violin_plot(sites[sites_to_plot], alpha_volumes[sites_to_plot], grd, tag="α", color="green", x_label="α")
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


        xsize = 300 + 70 * length(sites[sites_to_plot])
        plot!(size=(xsize, 300), margins=1.5Plots.cm, grid=false, legendfontsize=8)
        savefig(analysis_name * "_site_omega_means.pdf")

    end

    return df
end
export difFUBAR_tabulate

#Must return enough to re-calculate detections etc
export difFUBAR
"""
    difFUBAR(seqnames, seqs, treestring, tags, outpath; <keyword arguments>)

Takes a tagged phylogeny and an alignment as input and performs difFUBAR analysis.
Returns `df, results_tuple` where `df` is a DataFrame of the detected sites and `results_tuple` is a tuple of the partial calculations needed to re-run `difFUBAR_tabulate`.
Consistent with the docs of [`difFUBAR_tabulate`](@ref), `results_tuple` stores `(alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors)`.

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
    df = difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid; tag_colors=tag_colors, verbosity=verbosity, exports=exports)

    #Return df, (tuple of partial calculations needed to re-run tablulate)
    return df, (alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors)
end



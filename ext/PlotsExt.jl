module PlotsExt

using CodonMolecularEvolution
using Plots
using Measures
using MolecularEvolution
using Phylo
using DataFrames
using CSV
function gridplot(grid::CodonMolecularEvolution.FUBARgrid, results::CodonMolecularEvolution.BayesianFUBARResults; title="")
    θ = results.posterior_mean
    p = scatter(grid.alpha_ind_vec, grid.beta_ind_vec, zcolor=θ, c=:darktest, colorbar=false,
        markersize=sqrt(length(grid.alpha_ind_vec)) / 3.5, markershape=:square, markerstrokewidth=0.0, size=(350, 350),
        label=:none, xticks=(1:length(grid.grid_values), round.(grid.grid_values, digits=3)), xrotation=90,
        yticks=(1:length(grid.grid_values), round.(grid.grid_values, digits=3)), margin=6Plots.mm,
        xlabel="α", ylabel="β", title=title)
    plot(p, 1:length(grid.grid_values), 1:length(grid.grid_values), color="grey", style=:dash, label=:none)
    return p
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

function CodonMolecularEvolution.FUBAR_omega_plot(param_means, tag_colors, pos_thresh, detections, num_sites)
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
function CodonMolecularEvolution.plot_tagged_phylo_tree(tree, tag_colors, tags, analysis_name; strip_tags_from_name = CodonMolecularEvolution.generate_tag_stripper(tags))
    #TODO: update plots in docs
    phylo_tree = get_phylo_tree(tree)
    tagging = [tag_colors[CodonMolecularEvolution.model_ind(n, tags)] for n in nodenameiter(phylo_tree)]
    for node in nodeiter(phylo_tree)
        renamenode!(phylo_tree, node, strip_tags_from_name(node.name))
    end
    #Warnings regarding marker- and linecolor also appear in the Phylo.jl docs example
    #Note: sometimes long leafnames are truncated/not visible in the plot
    pl = plot(phylo_tree,
        showtips = true, tipfont = 6, markercolor = tagging, linecolor = tagging, markerstrokewidth = 0, size = (600, (120 + length(getleaflist(tree)) * 8)))
    savefig_tweakSVG(analysis_name * "_tagged_input_tree.svg", pl)
end
function CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid; tag_colors=DIFFUBAR_TAG_COLORS, verbosity=1, sites_to_plot=nothing, exports=true)
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

        site_counts_ω1 = CodonMolecularEvolution.collapse_counts(ω1, alloc_grid[:, site], cases=omegagrid)
        site_counts_ω2 = CodonMolecularEvolution.collapse_counts(ω2, alloc_grid[:, site], cases=omegagrid)
        site_counts_alphas = CodonMolecularEvolution.collapse_counts(alphas, alloc_grid[:, site], cases=alphagrid)

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
        CodonMolecularEvolution.FUBAR_omega_plot(param_means, tag_colors, pos_thresh, detections, num_sites)


        xsize = 300 + 70 * length(sites[sites_to_plot])
        plot!(size=(xsize, 300), margins=1.5Plots.cm, grid=false, legendfontsize=8)
        savefig(analysis_name * "_site_omega_means.pdf")

    end

    return df
end
function violin_plots(grid::CodonMolecularEvolution.FUBARgrid, results::CodonMolecularEvolution.BayesianFUBARResults;
    posterior_threshold=0.95,
    volume_scaling=1.0,
    plot_positive=true,
    plot_purifying=true)

    # Extract grid values for the x-axis
    grid_values = grid.alpha_vec[1:Int(sqrt(length(grid.alpha_vec)))]
    grd = round.(grid_values, digits=3)

    # Find sites with significant positive selection
    sites_positive = findall(results.positive_posteriors .> posterior_threshold)
    num_positive = length(sites_positive)
    println("$(num_positive) sites with positive selection above threshold.")

    p_positive = nothing
    p_purifying = nothing

    if num_positive > 0 && plot_positive
        p_positive = plot()
        # Calculate scaling factor for distributions
        s = 0.5 / max(maximum(results.posterior_alpha[:, sites_positive]),
            maximum(results.posterior_beta[:, sites_positive]))

        # Plot alpha distributions for positively selected sites
        FUBAR_violin_plot(sites_positive,
            [s .* volume_scaling .* results.posterior_alpha[:, [i]] for i in sites_positive],
            grd, tag="α", color="blue", legend_ncol=2, vertical_ind=nothing)

        # Plot beta distributions for positively selected sites
        FUBAR_violin_plot(sites_positive,
            [s .* volume_scaling .* results.posterior_beta[:, [i]] for i in sites_positive],
            grd, tag="β", color="red", legend_ncol=2, vertical_ind=nothing)

        plot!(size=(400, num_positive * 17 + 300), grid=false, margin=15Plots.mm)
    end

    # Find sites with significant purifying selection
    sites_purifying = findall(results.purifying_posteriors .> posterior_threshold)
    num_purifying = length(sites_purifying)
    println("$(num_purifying) sites with purifying selection above threshold.")

    if num_purifying > 0 && plot_purifying
        p_purifying = plot()
        # Calculate scaling factor for distributions
        s = 0.5 / max(maximum(results.posterior_alpha[:, sites_purifying]),
            maximum(results.posterior_beta[:, sites_purifying]))

        # Plot alpha distributions for purifying selected sites
        FUBAR_violin_plot(sites_purifying,
            [s .* volume_scaling .* results.posterior_alpha[:, [i]] for i in sites_purifying],
            grd, tag="α", color="blue", legend_ncol=2, vertical_ind=nothing)

        # Plot beta distributions for purifying selected sites
        FUBAR_violin_plot(sites_purifying,
            [s .* volume_scaling .* results.posterior_beta[:, [i]] for i in sites_purifying],
            grd, tag="β", color="red", legend_ncol=2, vertical_ind=nothing)

        plot!(size=(400, num_purifying * 17 + 300), grid=false, margin=15Plots.mm)
    end

    return p_positive, p_purifying
end


function CodonMolecularEvolution.plot_fubar_results(method::CodonMolecularEvolution.DefaultBayesianFUBARMethod, results::CodonMolecularEvolution.BayesianFUBARResults, grid::CodonMolecularEvolution.FUBARgrid; analysis_name="bayesian_analysis", write=false)
    posterior_mean_plot = gridplot(grid, results)
    positive_violin_plot, purifying_violin_plot = violin_plots(grid, results)
    if write
        savefig(posterior_mean_plot, analysis_name*"_posterior_mean.pdf")
        if !isnothing(positive_violin_plot)
            savefig(positive_violin_plot, analysis_name*"_positive_violin.pdf")
        end
        if !isnothing(purifying_violin_plot)
            savefig(purifying_violin_plot, analysis_name*"_purifying_violin.pdf")
        end
    end
    return posterior_mean_plot, violin_plots
end

function plot_skbdi_mixing(results::CodonMolecularEvolution.BayesianFUBARResults, grid::CodonMolecularEvolution.FUBARgrid, analysis_name)
    
    println(length(results.theta_chain[1]))

end

function CodonMolecularEvolution.plot_fubar_results(method::CodonMolecularEvolution.SKBDIFUBAR, results::CodonMolecularEvolution.BayesianFUBARResults, grid::CodonMolecularEvolution.FUBARgrid; analysis_name="skbdi_analysis", write=false, diagnostics = true)
    plot_skbdi_mixing(results, grid, analysis_name)
    CodonMolecularEvolution.plot_fubar_results(CodonMolecularEvolution.DefaultBayesianFUBARMethod(), results, grid, analysis_name=analysis_name, write=write)
end
function CodonMolecularEvolution.plot_fubar_results(method::CodonMolecularEvolution.DirichletFUBAR, results::CodonMolecularEvolution.BayesianFUBARResults, grid::CodonMolecularEvolution.FUBARgrid; analysis_name="dirichlet_analysis", write=false, diagnostics = false)
    CodonMolecularEvolution.plot_fubar_results(CodonMolecularEvolution.DefaultBayesianFUBARMethod(), results, grid, analysis_name=analysis_name, write=write)
end

end
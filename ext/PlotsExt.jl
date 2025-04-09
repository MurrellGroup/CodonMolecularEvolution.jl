module PlotsExt

using CodonMolecularEvolution
using Plots

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
    println("Yay we executed the loaded plotting fn")
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
function CodonMolecularEvolution.plot_fubar_results(method::CodonMolecularEvolution.SKBDIFUBAR, results::CodonMolecularEvolution.BayesianFUBARResults, grid::CodonMolecularEvolution.FUBARgrid; analysis_name="bayesian_analysis", write=false)
    print("Executing plotting skbdi")
    CodonMolecularEvolution.plot_fubar_results(CodonMolecularEvolution.DefaultBayesianFUBARMethod(), results, grid, analysis_name=analysis_name, write=write)
end
function CodonMolecularEvolution.plot_fubar_results(method::CodonMolecularEvolution.DirichletFUBAR, results::CodonMolecularEvolution.BayesianFUBARResults, grid::CodonMolecularEvolution.FUBARgrid; analysis_name="bayesian_analysis", write=false)
    println("Executing plotting dirichlet")
    CodonMolecularEvolution.plot_fubar_results(CodonMolecularEvolution.DefaultBayesianFUBARMethod(), results, grid, analysis_name=analysis_name, write=write)
end

end
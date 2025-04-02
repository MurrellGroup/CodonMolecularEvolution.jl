module PlotsExt

using CodonMolecularEvolution
using Plots

function CodonMolecularEvolution.gridplot(grid::CodonMolecularEvolution.FUBARgrid, θ; title="")
    scatter(grid.alpha_ind_vec, grid.beta_ind_vec, zcolor=θ, c=:darktest,
        markersize=sqrt(length(grid.alpha_ind_vec)) / 3.5, markershape=:square, markerstrokewidth=0.0, size=(400, 350),
        label=:none, xticks=(1:length(grid.grid_values), round.(grid.grid_values, digits=3)), xrotation=90,
        yticks=(1:length(grid.grid_values), round.(grid.grid_values, digits=3)), margin=6Plots.mm,
        xlabel="α", ylabel="β", title=title, colorbar=true, right_margin=12Plots.mm)
    plot!(1:length(grid.grid_values), 1:length(grid.grid_values), color="grey", style=:dash, label=:none)
end

function CodonMolecularEvolution.violin_plot_sites(grid::CodonMolecularEvolution.FUBARgrid, posteriors, posterior_alpha::Matrix,
    posterior_beta::Matrix;
    volume_scaling=1.0,
    posterior_threshold=0.95)

    grid_values = grid.beta_vec[1:Int(sqrt(length(grid.beta_vec)))]
    grd = round.(grid_values, digits=3)

    sites_to_plot = findall(posteriors .> posterior_threshold)
    num_plot = length(sites_to_plot)
    if num_plot == 0
        return nothing
    end
    violin_plots = plot()
    
    s = 0.5 / max(maximum(posterior_alpha[:, sites_to_plot]), maximum(posterior_beta[:, sites_to_plot]))
    FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_alpha[:, i:i] for i in sites_to_plot], grd, tag="α", color="blue", legend_ncol=2, vertical_ind=nothing)
    FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_beta[:, i:i] for i in sites_to_plot], grd, tag="β", color="red", legend_ncol=2, vertical_ind=nothing)
    plot!(size=(400, num_plot * 17 + 300), grid=false, margin=15Plots.mm)
    return violin_plots
end

function CodonMolecularEvolution.plot_FUBAR_posterior(grid::CodonMolecularEvolution.FUBARgrid, posterior::CodonMolecularEvolution.FUBARPosterior;
    posterior_threshold=0.95,
    volume_scaling=1.0)

    positive_violin_plots = violin_plot_sites(grid, posterior.positive_posteriors,
        posterior.posterior_alpha, posterior.posterior_beta,
        volume_scaling=volume_scaling,
        posterior_threshold=posterior_threshold)

    purifying_violin_plots = violin_plot_sites(grid, posterior.purifying_posteriors,
        posterior.posterior_alpha, posterior.posterior_beta,
        volume_scaling=volume_scaling,
        posterior_threshold=posterior_threshold)
        
    posterior_mean_plot = gridplot(grid, posterior.posterior_mean, title="Posterior mean θ")

    return positive_violin_plots, purifying_violin_plots, posterior_mean_plot
end

function CodonMolecularEvolution.plot_means_and_posteriors(param_means_vec, detections_vec)
    plots = []
    titles = String[]
    #Plot means
    categories = ["α", "ω1", "ω2"]
    for (i, category) in enumerate(categories)
        mean_codon = [p[i] for p in param_means_vec[1]]
        mean_nuc_codon = [p[i] for p in param_means_vec[2]]
        x = 0.0:0.01:maximum(mean_codon)
        p = plot(x, x, color=:black, label=false, aspect_ratio=:equal)
        scatter!(mean_codon, mean_nuc_codon, color=:blue, label=false, aspect_ratio=:equal)
        push!(titles, category * " means")
        push!(plots, p)
    end

    #Plot posterior probabilities
    categories = ["P(ω1 > ω2)", "P(ω2 > ω1)", "P(ω1 > 1)", "P(ω2 > 1)"]
    p2 = plot(layout = (1, length(categories)), thickness_scaling = 0.5)
    for (i, category) in enumerate(categories)
        codon = [p[i] for p in detections_vec[1]]
        nuc_codon = [p[i] for p in detections_vec[2]]
        x = 0.0:0.01:1.0
        p = plot(x, x, color=:black, label=false, aspect_ratio=:equal)
        scatter!(codon, nuc_codon, color=:red, label=false, aspect_ratio=:equal)
        push!(titles, category * " posterior")
        push!(plots, p)
    end
    l = @layout([° ° ° _; ° ° ° °])
    return plot(plots..., layout = l, size = (800, 400), title=reshape(titles, 1, 7), 
               plot_title="Codon vs. Nucleotide+Codon model fit", 
               xlabel="Codon", ylabel="Nucleotide+Codon", 
               left_margin=50px, bottom_margin=50px, thickness_scaling = 0.5)
end

function CodonMolecularEvolution.plot_trace(samples, param, LP, outpath; transform = x->x, title = "", ylabel = "")
    p = plot([transform(LP.unflatten(s.z.θ)[param]) for s in samples], 
         size = (700,350), 
         xlabel = "Iterations", 
         title = title, 
         ylabel = ylabel, 
         legend = :none)
    save_plot(p, outpath)
    return p
end

function CodonMolecularEvolution.core_plots(ℓπ, samples, posterior_mean_θ, f, analysis_name, n_adapts; plots = true)
    if plots
        anim = @animate for i ∈ (n_adapts+1):1:length(samples)
            gridplot(f.alpha_ind_vec,f.beta_ind_vec,f.grid_values,thetas(ℓπ, samples[i].z.θ))
        end
        gif(anim, analysis_name * "_posterior_θ_samples.mp4", fps = 15)
        
        log_posterior_plot = plot(f.LL_offset .+ [samples[i].stat.log_density for i in (n_adapts+1):length(samples)], 
             label = "Log Posterior", 
             size = (700,350), 
             xlabel = "Iterations")
        save_plot(log_posterior_plot, analysis_name * "_log_posterior_trace.pdf")
        
        log_posterior_burnin_plot = plot(f.LL_offset .+ [samples[i].stat.log_density for i in 1:length(samples)], 
             label = "Log Posterior", 
             size = (700,350), 
             xlabel = "Iterations")
        save_plot(log_posterior_burnin_plot, analysis_name * "_log_posterior_trace_with_burnin.pdf")
        
        all_theta_chains = stack([thetas(ℓπ, samples[s].z.θ) for s in (n_adapts+1):length(samples)])
        inds_to_plot = sortperm(posterior_mean_θ, rev = true)[(1:(Int(sqrt(length(f.alpha_ind_vec))))).^2]
        theta_trace_plot = plot(all_theta_chains[inds_to_plot,:]',
             legend = :none, 
             alpha = 0.5, 
             size = (700,350), 
             xlabel = "Iterations")
        save_plot(theta_trace_plot, analysis_name * "_posterior_θ_trace.pdf")
    end
end

function CodonMolecularEvolution.other_plots(LP::LogPosteriorRFF{Float64, M}, samples, f, analysis_name, n_adapts, posterior_mean_θ; plots = true, verbosity = 1) where M<:WeightedHMC_FUBAR
    prior_possel = 1 - cdf(LP.priors.φpre,0)
    poschain = [LP.unflatten(s.z.θ).φpre > 0 for s in samples[n_adapts+1:end]]
    if all(poschain)
        pospos = 1.0
        bayesfactor = Inf
    else
        pospos = sum(poschain/length(poschain))
        bayesfactor = pospos/(1-pospos) * (1-prior_possel)/prior_possel
    end
    ppstr = "P(φ>0|data)=$(round(pospos,sigdigits=4)). BF = $(round(bayesfactor,sigdigits=4))"
    
    if plots
        φpre_trace = plot_trace(samples[n_adapts+1:end], :φpre, LP, nothing; 
                  transform = x->x, 
                  title = ppstr, 
                  ylabel = "φpre")
        save_plot(φpre_trace, analysis_name*"_φpre_trace.pdf")
        
        φ_trace = plot_trace(samples[n_adapts+1:end], :φpre, LP, nothing; 
                  transform = transf(M()), 
                  title = ppstr, 
                  ylabel = "φ")
        save_plot(φ_trace, analysis_name*"_φ_trace.pdf")
    end
    
    smoothFUBAR_mixing = calculate_smoothFUBAR_mixing(LP, samples, n_adapts)
    return (θ = posterior_mean_θ, φ_crossings = smoothFUBAR_mixing, global_posterior_probability = pospos, global_bayes_factor = bayesfactor)
end

function save_plot(p, filename)
    if !isnothing(filename)
        savefig(p, filename)
    end
    return p
end

# Override the default implementation when Plots is available
function CodonMolecularEvolution.maybe_save_plots(plots_to_save::NamedTuple, save::Bool)
    if !save
        return nothing
    end
    
    for (plot, filename) in values(plots_to_save)
        if !isnothing(plot)
            save_plot(plot, filename)
        end
    end
end

end 
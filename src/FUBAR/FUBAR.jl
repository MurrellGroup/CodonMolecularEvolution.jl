#Need to roll this into original FUBAR
struct FUBARgrid{T}
    grid_values::Vector{T}
    alpha_vec::Vector{T}
    beta_vec::Vector{T}
    alpha_ind_vec::Vector{Int64}
    beta_ind_vec::Vector{Int64}
    cond_lik_matrix::Matrix{T}
    LL_offset::T
    sites::Int64
end


function gridplot(alpha_ind_vec,beta_ind_vec,grid_values,θ; title = "")
    scatter(alpha_ind_vec,beta_ind_vec, zcolor = θ, c = :darktest, colorbar = false,
    markersize = sqrt(length(alpha_ind_vec))/3.5, markershape=:square, markerstrokewidth=0.0, size=(350,350),
    label = :none, xticks = (1:length(grid_values), round.(grid_values,digits = 3)), xrotation = 90,
    yticks = (1:length(grid_values), round.(grid_values,digits = 3)), margin=6Plots.mm,
    xlabel = "α", ylabel = "β", title = title)
    plot!(1:length(grid_values),1:length(grid_values),color = "grey", style = :dash, label = :none)
end

function FUBAR_init(outpath_and_file_prefix, treestring; verbosity=1, exports=true, disable_binarize=false, ladderize_tree = false)
    analysis_name = outpath_and_file_prefix
    splt = splitpath(analysis_name)[1:end-1]
    if length(splt) > 0
        exports && mkpath(joinpath(splt))
    end
    tree = gettreefromnewick(treestring, FelNode, disable_binarize=disable_binarize)
    if ladderize_tree
        MolecularEvolution.ladderize!(tree)
    end
    verbosity > 0 && println("Step 1: Initialization.")
    return tree, analysis_name
end

sites(p::LazyPartition{CodonPartition}) = p.memoryblocks[1].sites
sites(p::CodonPartition) = p.sites


function FUBAR_grid(tree, GTRmat, F3x4_freqs, code; verbosity=1)
    verbosity > 0 && println("Step 3: Calculating conditional likelihoods.")

    grid_values = 10 .^ (-1.35:0.152:1.6) .- 0.0423174293933042 #Possibly should be an argument

    LL_matrix = zeros(length(grid_values)^2,sites(tree.message[1]));
    alpha_vec = zeros(length(grid_values)^2);
    alpha_ind_vec = zeros(Int64,length(grid_values)^2);
    beta_vec = zeros(length(grid_values)^2);
    beta_ind_vec = zeros(Int64,length(grid_values)^2);

    i = 1
    for (a,alpha) in enumerate(grid_values)
        for (b,beta) in enumerate(grid_values)
            alpha_vec[i],beta_vec[i] = alpha, beta
            alpha_ind_vec[i], beta_ind_vec[i] = a,b
            m = DiagonalizedCTMC(MolecularEvolution.MG94_F3x4(alpha, beta, GTRmat, F3x4_freqs))
            felsenstein!(tree,m)
            combine!(tree.message[1],tree.parent_message[1])
            LL_matrix[i,:] .= MolecularEvolution.site_LLs(tree.message[1])
            i += 1
        end
    end
    maxi_shift = maximum(LL_matrix,dims = 1)
    prob_matrix = exp.(LL_matrix .- maxi_shift)
    sum_shift = sum(prob_matrix,dims = 1)
    prob_matrix ./= sum_shift
    LLO = sum(maxi_shift .+ log.(sum_shift))
    return FUBARgrid(grid_values, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec, prob_matrix, LLO, size(prob_matrix,2))
end

function FUBAR_fitEM(con_lik_matrix, iters, conc; verbosity=1)
    verbosity > 0 && println("Step 4: Model fitting.")
    L = size(con_lik_matrix,1)
    LDAθ = weightEM(con_lik_matrix, ones(L)./L, conc = conc, iters = iters);
    return LDAθ
end

function FUBAR_tabulate_from_θ(θ, f::FUBARgrid, analysis_name; posterior_threshold = 0.95, volume_scaling = 1.0, verbosity = 1)
    verbosity > 0 && println("Step 5: Tabulating results and saving plots.")
    con_lik_matrix, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec = f.cond_lik_matrix, f.alpha_vec, f.beta_vec, f.alpha_ind_vec, f.beta_ind_vec
    pos_filt = beta_ind_vec .> alpha_ind_vec
    pur_filt = beta_ind_vec .< alpha_ind_vec
    weighted_mat = con_lik_matrix .* θ
    weighted_mat ./= sum(weighted_mat,dims = 1)
    positive_posteriors = sum(weighted_mat[pos_filt,:],dims = 1)[:]
    purifying_posteriors = sum(weighted_mat[pur_filt,:],dims = 1)[:]
    beta_pos_mean = sum(weighted_mat .* beta_vec, dims = 1)[:]
    alpha_pos_mean = sum(weighted_mat .* alpha_vec, dims = 1)[:]

    weighted_sites = reshape(weighted_mat, 20,20,:);
    posterior_alpha = sum(weighted_sites, dims = 1)[1,:,:]
    posterior_beta = sum(weighted_sites, dims = 2)[:,1,:]

    grid_values = beta_vec[1:Int(sqrt(length(beta_vec)))]
    grd = round.(grid_values, digits = 3)

    sites_to_plot = findall(positive_posteriors .> posterior_threshold)
   #= num_plot = length(sites_to_plot)
    if num_plot > 0
        verbosity > 0 && println("$num_plot sites with positive selection above threshold.")
        plot()
        s = 0.5/max(maximum(posterior_alpha[:,sites_to_plot]),maximum(posterior_beta[:,sites_to_plot]))
        FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_alpha[:,[i]] for i in sites_to_plot], grd, tag="α", color="blue", legend_ncol=2, vertical_ind = nothing)
        FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_beta[:,[i]] for i in sites_to_plot], grd, tag="β", color="red", legend_ncol=2, vertical_ind = nothing)
        plot!(size=(400, num_plot * 17 + 300), grid=false, margin=15Plots.mm)
        savefig(analysis_name * "_violin_positive.pdf")
    else
        verbosity > 0 && println("No sites with positive selection above threshold.")
    end

    sites_to_plot = findall(purifying_posteriors .> posterior_threshold)
    num_plot = length(sites_to_plot)
    if num_plot > 0
        verbosity > 0 && println("$num_plot sites with purifying selection above threshold.")
        plot()
        s = 0.5/max(maximum(posterior_alpha[:,sites_to_plot]),maximum(posterior_beta[:,sites_to_plot]))
        FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_alpha[:,[i]] for i in sites_to_plot], grd, tag="α", color="blue", legend_ncol=2, vertical_ind = nothing)
        FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_beta[:,[i]] for i in sites_to_plot], grd, tag="β", color="red", legend_ncol=2, vertical_ind = nothing)        
        plot!(size=(400, num_plot * 17 + 300), grid=false, margin=15Plots.mm)
        savefig(analysis_name * "_violin_purifying.pdf")
    else
        verbosity > 0 && println("No sites with purifying selection above threshold.")
    end

    gridplot(alpha_ind_vec,beta_ind_vec,grid_values,θ; title = "Posterior mean θ")
    savefig(analysis_name * "_θ.pdf")
    =#
    df_results = DataFrame(
        site = 1:size(con_lik_matrix,2),
        positive_posterior = positive_posteriors,
        purifying_posterior = purifying_posteriors,
        beta_pos_mean = beta_pos_mean,
        alpha_pos_mean = alpha_pos_mean,
    )
    CSV.write(analysis_name * "_results.csv", df_results)
    return df_results
end

function FUBAR_plot_from_θ(θ, f::FUBARgrid, analysis_name; posterior_threshold = 0.95, volume_scaling = 1.0, verbosity = 1)
    verbosity > 0 && println("Step 5: Tabulating results and saving plots.")
    con_lik_matrix, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec = f.cond_lik_matrix, f.alpha_vec, f.beta_vec, f.alpha_ind_vec, f.beta_ind_vec
    pos_filt = beta_ind_vec .> alpha_ind_vec
    pur_filt = beta_ind_vec .< alpha_ind_vec
    weighted_mat = con_lik_matrix .* θ
    weighted_mat ./= sum(weighted_mat,dims = 1)
    positive_posteriors = sum(weighted_mat[pos_filt,:],dims = 1)[:]
    purifying_posteriors = sum(weighted_mat[pur_filt,:],dims = 1)[:]
    beta_pos_mean = sum(weighted_mat .* beta_vec, dims = 1)[:]
    alpha_pos_mean = sum(weighted_mat .* alpha_vec, dims = 1)[:]

    weighted_sites = reshape(weighted_mat, 20,20,:);
    posterior_alpha = sum(weighted_sites, dims = 1)[1,:,:]
    posterior_beta = sum(weighted_sites, dims = 2)[:,1,:]

    grid_values = beta_vec[1:Int(sqrt(length(beta_vec)))]
    grd = round.(grid_values, digits = 3)

    sites_to_plot = findall(positive_posteriors .> posterior_threshold)
    num_plot = length(sites_to_plot)
    if num_plot > 0
        verbosity > 0 && println("$num_plot sites with positive selection above threshold.")
        plot()
        s = 0.5/max(maximum(posterior_alpha[:,sites_to_plot]),maximum(posterior_beta[:,sites_to_plot]))
        FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_alpha[:,[i]] for i in sites_to_plot], grd, tag="α", color="blue", legend_ncol=2, vertical_ind = nothing)
        FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_beta[:,[i]] for i in sites_to_plot], grd, tag="β", color="red", legend_ncol=2, vertical_ind = nothing)
        plot!(size=(400, num_plot * 17 + 300), grid=false, margin=15Plots.mm)
        savefig(analysis_name * "_violin_positive.pdf")
    else
        verbosity > 0 && println("No sites with positive selection above threshold.")
    end

    sites_to_plot = findall(purifying_posteriors .> posterior_threshold)
    num_plot = length(sites_to_plot)
    if num_plot > 0
        verbosity > 0 && println("$num_plot sites with purifying selection above threshold.")
        plot()
        s = 0.5/max(maximum(posterior_alpha[:,sites_to_plot]),maximum(posterior_beta[:,sites_to_plot]))
        FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_alpha[:,[i]] for i in sites_to_plot], grd, tag="α", color="blue", legend_ncol=2, vertical_ind = nothing)
        FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_beta[:,[i]] for i in sites_to_plot], grd, tag="β", color="red", legend_ncol=2, vertical_ind = nothing)        
        plot!(size=(400, num_plot * 17 + 300), grid=false, margin=15Plots.mm)
        savefig(analysis_name * "_violin_purifying.pdf")
    else
        verbosity > 0 && println("No sites with purifying selection above threshold.")
    end

    gridplot(alpha_ind_vec,beta_ind_vec,grid_values,θ; title = "Posterior mean θ")
    savefig(analysis_name * "_θ.pdf")

    
end


#Packaging "everything before the conditional likelihoods"
function FUBAR_init2grid(seqnames, seqs, treestring, outpath;
    pos_thresh=0.95, verbosity=1, exports=true, code=MolecularEvolution.universal_code, optimize_branch_lengths=false)
    analysis_name = outpath
    tree, analysis_name = FUBAR_init(analysis_name, treestring, exports=exports, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit_2steps(seqnames, seqs, tree, x -> x, code, verbosity=verbosity, optimize_branch_lengths=optimize_branch_lengths)
    return FUBAR_grid(tree, GTRmat, F3x4_freqs, code, verbosity=verbosity)
end

export FUBAR_init2grid, FUBAR_tabulate_from_θ

function FUBAR(seqnames, seqs, treestring, outpath;
    pos_thresh=0.95, verbosity=1, exports=true, code=MolecularEvolution.universal_code, optimize_branch_lengths=false,
    method = (sampler = :DirichletEM, concentration = 0.5, iterations = 2500))
    fubar = FUBAR_init2grid(seqnames, seqs, treestring, outpath,
    pos_thresh=pos_thresh, verbosity=verbosity, exports=exports, code=code, optimize_branch_lengths=optimize_branch_lengths)
    
    con_lik_matrix, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec, LL_offset = fubar.cond_lik_matrix, fubar.alpha_vec, fubar.beta_vec, fubar.alpha_ind_vec, fubar.beta_ind_vec, fubar.cond_lik_matrix, fubar.LL_offset

    #if method.sampler == :DirichletEM
        θ = FUBAR_fitEM(con_lik_matrix, method.iterations, method.concentration)
    #end
    
    df_results = FUBAR_tabulate_from_θ(θ, fubar, outpath, posterior_threshold = pos_thresh, verbosity = verbosity)
    #Return df, (tuple of partial calculations needed to re-run tablulate)
    return df_results, (θ, fubar)
end

function FUBAR_precomputed_f(f, outpath;
    pos_thresh=0.95, verbosity=1, exports=true, code=MolecularEvolution.universal_code, optimize_branch_lengths=false,
    method = (sampler = :DirichletEM, concentration = 0.5, iterations = 2500))
   
    
    con_lik_matrix, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec, LL_offset = f.cond_lik_matrix, f.alpha_vec, f.beta_vec, f.alpha_ind_vec, f.beta_ind_vec, f.cond_lik_matrix, f.LL_offset

    #if method.sampler == :DirichletEM
        θ = FUBAR_fitEM(con_lik_matrix, method.iterations, method.concentration)
    #end
    
    df_results = FUBAR_tabulate_from_θ(θ, f, outpath, posterior_threshold = pos_thresh, verbosity = verbosity)
    #Return df, (tuple of partial calculations needed to re-run tablulate)
    return df_results, (θ, f)
end

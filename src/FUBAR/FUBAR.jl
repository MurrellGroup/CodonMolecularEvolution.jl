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
    grid_function::Function
    LL_matrix::Matrix{T}
end


function gridplot(alpha_ind_vec,beta_ind_vec,grid_values,θ; title = "")
    scatter(alpha_ind_vec,beta_ind_vec, zcolor = θ, c = :darktest, colorbar = false,
    markersize = sqrt(length(alpha_ind_vec))/3.5, markershape=:square, markerstrokewidth=0.0, size=(350,350),
    label = :none, xticks = (1:length(grid_values), round.(grid_values,digits = 3)), xrotation = 90,
    yticks = (1:length(grid_values), round.(grid_values,digits = 3)), margin=6Plots.mm,
    xlabel = "α", ylabel = "β", title = title)
    plot!(1:length(grid_values),1:length(grid_values),color = "grey", style = :dash, label = :none)
end

function init_path(analysis_name)
    splt = splitpath(analysis_name)[1:end-1]
    if length(splt) > 0
        mkpath(joinpath(splt))
    end
end

function FUBAR_init(treestring; verbosity=1, exports=true, disable_binarize=false, ladderize_tree = false)

    tree = gettreefromnewick(treestring, FelNode, disable_binarize=disable_binarize)
    if ladderize_tree
        MolecularEvolution.ladderize!(tree)
    end
    verbosity > 0 && println("Step 1: Initialization.")
    return tree
end

sites(p::LazyPartition{CodonPartition}) = p.memoryblocks[1].sites
sites(p::CodonPartition) = p.sites


function FUBAR_grid(tree, GTRmat, F3x4_freqs, code; verbosity=1, grid_function = x -> 10^(x/6.578947368421053 - 1.502) - 0.0423174293933042, num_grid_points = 20)
    verbosity > 0 && println("Step 3: Calculating conditional likelihoods.")

    grid_values = grid_function.(1:num_grid_points)

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
    return FUBARgrid(grid_values, alpha_vec, beta_vec, alpha_ind_vec, beta_ind_vec, prob_matrix, LLO, size(prob_matrix,2), grid_function, LL_matrix)
end

function FUBAR_fitEM(con_lik_matrix, iters, conc; verbosity=1)
    verbosity > 0 && println("Step 4: Model fitting.")
    L = size(con_lik_matrix,1)
    LDAθ = weightEM(con_lik_matrix, ones(L)./L, conc = conc, iters = iters);
    return LDAθ
end

function FUBAR_tabulate_from_θ(θ, f::FUBARgrid, analysis_name; posterior_threshold = 0.95, volume_scaling = 1.0, verbosity = 1, plots = true)
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
    #Note: summing over first dim is summing over beta
    posterior_alpha = sum(weighted_sites, dims = 1)[1,:,:]
    posterior_beta = sum(weighted_sites, dims = 2)[:,1,:]

    grid_values = beta_vec[1:Int(sqrt(length(beta_vec)))]
    grd = round.(grid_values, digits = 3)

    sites_to_plot = findall(positive_posteriors .> posterior_threshold)
    num_plot = length(sites_to_plot)
    verbosity > 0 && println("$(num_plot) sites with positive selection above threshold.")
    if (num_plot > 0) && plots
        plot()
        s = 0.5/max(maximum(posterior_alpha[:,sites_to_plot]),maximum(posterior_beta[:,sites_to_plot]))
        FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_alpha[:,[i]] for i in sites_to_plot], grd, tag="α", color="blue", legend_ncol=2, vertical_ind = nothing)
        FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_beta[:,[i]] for i in sites_to_plot], grd, tag="β", color="red", legend_ncol=2, vertical_ind = nothing)
        plot!(size=(400, num_plot * 17 + 300), grid=false, margin=15Plots.mm)
        savefig(analysis_name * "_violin_positive.pdf")
    end
    
    sites_to_plot = findall(purifying_posteriors .> posterior_threshold)
    num_plot = length(sites_to_plot)
    verbosity > 0 && println("$(num_plot) sites with purifying selection above threshold.")
    if (num_plot > 0) && plots
        plot()
        s = 0.5/max(maximum(posterior_alpha[:,sites_to_plot]),maximum(posterior_beta[:,sites_to_plot]))
        FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_alpha[:,[i]] for i in sites_to_plot], grd, tag="α", color="blue", legend_ncol=2, vertical_ind = nothing)
        FUBAR_violin_plot(sites_to_plot, [s .* volume_scaling .* posterior_beta[:,[i]] for i in sites_to_plot], grd, tag="β", color="red", legend_ncol=2, vertical_ind = nothing)        
        plot!(size=(400, num_plot * 17 + 300), grid=false, margin=15Plots.mm)
        savefig(analysis_name * "_violin_purifying.pdf")
    end
    
    if plots
        gridplot(alpha_ind_vec,beta_ind_vec,grid_values,θ; title = "Posterior mean θ")
        savefig(analysis_name * "_θ.pdf")
    end
    
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

#Packaging "everything before the conditional likelihoods"
function alphabetagrid(seqnames::Vector{String}, seqs, treestring::String;
    verbosity=1, code=MolecularEvolution.universal_code, optimize_branch_lengths=false)
    tree = FUBAR_init(treestring, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit_2steps(seqnames, seqs, tree, x -> x, code, verbosity=verbosity, optimize_branch_lengths=optimize_branch_lengths)
    return FUBAR_grid(tree, GTRmat, F3x4_freqs, code, verbosity=verbosity)
end

export alphabetagrid, FUBAR_tabulate_from_θ

#Need to make this match the smoothFUBAR setup with the first argument controlling the method (EM, Gibbs, etc)
function FUBAR(f::FUBARgrid, outpath;
    pos_thresh=0.95, verbosity=1, exports=true, code=MolecularEvolution.universal_code, optimize_branch_lengths=false,
    method = (sampler = :DirichletEM, concentration = 0.5, iterations = 2500), plots = true)
    exports && init_path(outpath)
    θ = FUBAR_fitEM(f.cond_lik_matrix, method.iterations, method.concentration, verbosity=verbosity)
    df_results = FUBAR_tabulate_from_θ(θ, f, outpath, posterior_threshold = pos_thresh, verbosity = verbosity, plots = plots)
    return df_results, (θ = θ, )
end

function alternating_maximize(f, a_bounds, b_bounds; final_tol = 1e-20, max_iters = 50)
    a = sum(a_bounds)/2
    b = sum(b_bounds)/2
    m = -f(a,b)
    for i in 1:3
        a = brents_method_minimize(x -> -f(x,b), a_bounds[1], a_bounds[2], identity, 1e-5)
        b = brents_method_minimize(x -> -f(a,x), b_bounds[1], b_bounds[2], identity, 1e-5)
    end
    m = -f(a,b)
    next_m = -Inf
    for i in 1:max_iters
        a = brents_method_minimize(x -> -f(x,b), a_bounds[1], a_bounds[2], identity, final_tol)
        b = brents_method_minimize(x -> -f(a,x), b_bounds[1], b_bounds[2], identity, final_tol)
        next_m = -f(a,b)
        if abs(next_m - m) < final_tol
            break;
        end
        m = next_m
    end
    return a,b,-m
end

#Generalize this to work with any grid dimensions!
function interpolating_LRS(grid)
    itp = interpolate(grid, BSpline(Cubic(Line(OnGrid()))))
    
    #Null model:
    ab = brents_method_minimize(x -> -itp(x,x), 1, 20, identity, 1e-20)
    LL_null = itp(ab,ab)

    #Alt model:
    a,b,LL_alt = alternating_maximize(itp, (1.0,20.0), (1.0,20.0))
    LRS = 2(LL_alt-LL_null)
    p = 1-cdf(Chisq(1), LRS)
    return p, (p_value = p, LRS = LRS, LL_alt = LL_alt, α_alt = a, β_alt = b, LL_null = LL_null, αβ_null = ab, sel = ifelse(p<0.05,ifelse(b>a, "Positive", "Purifying"), ""))
end

#Frequentist interpolating fixed-effects FUBAR
#This needs to have exports like regular FUBAR - plots etc
#We could plot LL surfaces for each site under selection, like this:
#=
    itp = interpolate((LLmatrix[:,:,190]'), BSpline(Cubic(Line(OnGrid()))));
    x = range(1, 20, length=200)
    y = range(1, 20, length=200)
    z = @. itp(x', y)
    contour(x, y, z, levels=30, color=:turbo, fill=true, linewidth = 0, colorbar = false, size = (400,400))
=#
function FIFE(f::FUBARgrid, outpath; verbosity=1, exports=true)
    exports && init_path(outpath)
    LLmatrix = reshape(f.LL_matrix, length(f.grid_values),length(f.grid_values),:) 
    #Note: dim1 = beta, dim2 = alpha, so we transpose going in:
    stats = [interpolating_LRS(LLmatrix[:,:,i]') for i in 1:size(LLmatrix, 3)]
    df_results = DataFrame([s[2] for s in stats])
    df_results.site = 1:size(LLmatrix, 3)
    df_results.α_alt .= f.grid_function.(df_results.α_alt)
    df_results.β_alt .= f.grid_function.(df_results.β_alt)
    df_results.αβ_null .= f.grid_function.(df_results.αβ_null)
    CSV.write(outpath * "_results.csv", df_results)
    return df_results
end

export FIFE

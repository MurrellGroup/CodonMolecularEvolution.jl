#A call should look like
#f = FLAVORgrid(seqnames, seqs, treestring)
#FLAVOR(f, outdir)

const FLAVOR_CAP = 1.0

struct FLAVORgrid{T}
    tr::Function
    trinv::Function
    mugrid::Vector{T}
    shapegrid::Vector{T}
    alphagrid::Vector{T}
    gridpoints::Vector{Tuple{T,T,T}}
    grid_dims::Vector{Int64}
    prob_matrix::Matrix{T}
    site_scalers::Vector{T}
end


function FLAVORgrid(seqnames::Vector{String}, seqs, treestring::String;
    verbosity=1, code=MolecularEvolution.universal_code, optimize_branch_lengths=false)
    tree = FUBAR_init(treestring, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit_2steps(seqnames, seqs, tree, x -> x, code, verbosity=verbosity, optimize_branch_lengths=optimize_branch_lengths)
    
    verbosity > 0 && println("Step 3: Calculating conditional likelihood grid.")

    tr(x) = max(0.0,10^x-0.05)
    trinv(x) =  log10(x+0.05)

    mugrid = gridsetup(0.01, 16.0, 8, trinv, tr) 
    shapegrid = gridsetup(0.05, 20, 6, trinv, tr)
    alphagrid = gridsetup(0.01, 10, 8, trinv, tr)

    #Tiny grid for testing:
    #mugrid = gridsetup(0.01, 16.0, 3, trinv, tr) 
    #shapegrid = gridsetup(0.05, 20, 3, trinv, tr)
    #alphagrid = gridsetup(0.01, 10, 3, trinv, tr)

    #mu, shape alpha
    gridpoints = [(mu,shape,alpha) for mu in mugrid for shape in shapegrid for alpha in alphagrid];
    grid_dims = length.([mugrid,shapegrid,alphagrid])

    LL_matrix = gamma_BWMM_grid(mugrid,shapegrid,alphagrid,grid_dims, GTRmat, F3x4_freqs, tree, verbosity = verbosity)

    maxi_shift = maximum(LL_matrix,dims = 1)
    prob_matrix = exp.(LL_matrix .- maxi_shift)
    sum_shift = sum(prob_matrix,dims = 1)
    prob_matrix ./= sum_shift
    site_scalers = maxi_shift[:] .+ log.(sum_shift[:])

    return FLAVORgrid(tr, trinv, mugrid, shapegrid, alphagrid, gridpoints, grid_dims, prob_matrix, site_scalers)
end

export FLAVORgrid


#Converts from Gamma parameters to discrete gamma approximation, using quantiles.
function gamma_slices(mu,shape,slices)
    c = (1/slices)/2
    return quantile(Gamma(shape,mu/shape),c:2*c:1-c)
end

#This will construct a time sequence and a P sequence that you need for the BWMM
#Note: When doing grid work, other alphas should be obtained without re-running this, by calling "rescale!(model,factor)"
#I recommend you structure the grid loop so that alpha is on the inner most loop
#We should write a cleaner version of this that lets you use weights as well
#Make a version that will work for arbitrary genetics codes
#Sets up the codon InterpolatedDiscreteModel, based on MG94F3x4
function omega_BWMM_matrix_sequence(alpha::Float64,omega_vec::Vector{Float64}, nuc_mat::Array{Float64,2}, F3x4::Array{Float64,2}; t = 0.001,n = 50, cap = n-15)
    eq_freqs = MolecularEvolution.F3x4_eq_freqs(F3x4)
    Ps = zeros(61,61,n)
    for o in omega_vec
        Ps .+= (MolecularEvolution.matrix_sequence(MolecularEvolution.MG94_F3x4(alpha,alpha*o,nuc_mat,F3x4),t,n,cap = cap) ./ length(omega_vec))
    end
    #Some special handling to force eq freqs for numerical infinity - unclear if this matters
    Ps[:,:,end] .= reshape(repeat(eq_freqs, inner=61, outer=1),61,61)
    ts = MolecularEvolution.t_sequence(t,n,cap = cap)
    return ts, Ps
end

#Model wrapper - returns models function that gets used in eg. felsenstein_up!
function construct_GammaBWMM(alpha, mean, shape, nuc_matrix, F3x4; num_parts = 20, capped=false)
    omega_vec = gamma_slices(mean,shape,num_parts)
    if capped
        clamp!(omega_vec,0.0,1.0)
    end
    ts,ps = omega_BWMM_matrix_sequence(alpha,omega_vec,nuc_matrix, F3x4)
    m = MolecularEvolution.InterpolatedDiscreteModel(ps, ts);
    return m
end


function gamma_BWMM_grid(mugrid,shapegrid,alphagrid,grid_dims, GTRmat, F3x4_freqs, tree::FelNode; verbosity = 1)
    #Calculates the conditional likelihood grid.
    num_sites = CodonMolecularEvolution.sites(tree.message[1])
    l = prod(grid_dims)*2
    log_con_lik_matrix = zeros(l,num_sites)
    row_ind = 1
    
    for capped in [false,true]
        verbosity > 0 && print("Capped: ", capped, " ")
        for mu in mugrid
            for shape in shapegrid
                #Construct a model with alpha=1.0
                m = construct_GammaBWMM(1.0, mu, shape, GTRmat, F3x4_freqs, capped=capped)
                #Copy starting "ts" from the model
                alpha1tvec = copy(m.tvec)
                for alpha in alphagrid
                    #Generate all alpha models by just shifting the interpolation grid points
                    #which is equivalent to scaling things by a constant
                    m.tvec= alpha1tvec ./ alpha
                    felsenstein!(tree,m)
                    combine!(tree.message[1],tree.parent_message[1])
                    log_con_lik_matrix[row_ind,:] .= MolecularEvolution.site_LLs(tree.message[1])
                    row_ind += 1
                end
            end
            verbosity > 0 && print(round(mu, sigdigits = 3), " ")
        end
    end
    return log_con_lik_matrix
end

function bayes_factor(posterior,prior)
    return (posterior/(1-posterior))*((1-prior)/prior)
end

function unflatten(θ,grid_dims)
    mubyshapebyalpha = zeros(grid_dims...);
    c = 1
    for i in 1:grid_dims[1]
        for j in 1:grid_dims[2]
            for k in 1:grid_dims[3]
                mubyshapebyalpha[i,j,k] = θ[c]
                c += 1
            end
        end
    end
    return mubyshapebyalpha
end


#Visualizing the omega distributions (as CDFs) for the specified grid.
#Note the x axis has been scaled with x^4
num2hex(n) = lpad(string(min(max(Int64(round(n)),0),255), base=16),2,'0')
RGBhex(r,g,b) = "#"*num2hex(r)*num2hex(g)*num2hex(b)

function plot_gamma!(mu,shape; up = 150.0, alpha = 0.1, color = "black", slices = 20)
    slice_pts = gamma_slices(mu,shape,slices)
    slice_ys = cdf.(Gamma(shape,mu/shape),slice_pts) 
    xs = (trinv(0.0):0.001:(trinv(up)))
    v = cdf.(Gamma(shape,mu/shape),tr.(xs))
    step = (trinv(1.0) .- trinv(0.0))/10
    plot!(xs,v,alpha = alpha, color = color, xticks = ([trinv(0):step:trinv(up);], round.(tr.([trinv(0):step:trinv(up);]),sigdigits = 4)), xrotation = 90)
    scatter!(trinv.(slice_pts),slice_ys,alpha = alpha, color = color, markerstrokewidth = 0, markersize = 2.0)
end

function plot_gamma_distributions(mugrid, shapegrid)
    pl = plot(legend = false, size = (700,300), margins = 10Plots.mm)
    for mu in mugrid
        for shape in shapegrid
            gammastd = std(Gamma(shape,mu/shape))
            plot_gamma!(mu,shape,color = RGBhex(gammastd*40,0,256-gammastd*40), alpha = min(gammastd+0.2,0.5))
        end
    end
    return pl
end

#Calculates the proportion of omega categories > 1.
function prop_pos(gp, capped)
    mu,shape,alpha = gp
    s = gamma_slices(mu,shape,20)
    if capped
        return 0.0
    end
    return mean(s .> 1.0)
end

#Get pdf of a Distribution that's been right-capped-off at cap
function cappedpdf(G::Distribution, x::AbstractVector{Float64}, cap::Float64, width::Float64; capped = true, trinv = identity)
    if !capped
        return pdf(G, x)
    end
    cap_ind = searchsortedlast(x, cap)

    y = zeros(length(x))
    @views y[1:cap_ind] .+= pdf(G, x[1:cap_ind])
    distribute_point_mass!(y, G, x, true, cap, width, trinv=trinv)
    return y
end

function distribute_point_mass!(y::AbstractVector{Float64}, G::Distribution, x::AbstractVector{Float64}, to_right::Bool, loc::Float64, width::Float64; trinv = identity)
    L = length(x)
    @assert L == length(y)
    indfind = to_right ? max(1, searchsortedlast(x, loc)) : min(L, searchsortedfirst(x, loc))
    left_mass = cdf(G, x[indfind])
    mass = to_right ? 1 - left_mass : left_mass

    indstop = indfind
    trinv_xfind = trinv(x[indfind])
    while  abs(trinv(x[indstop]) - trinv_xfind) < width
        indstop = to_right ? max(1, indstop - 1) : min(L, indstop + 1)
        (indstop == 1 || indstop == L) && break
    end
    step = to_right ? -1 : 1
    @views y[indfind:step:indstop] .+= mass / abs(x[indstop] - x[indfind])
end

#Should only do this for detecs, because it's expensive
#DRY, see FUBAR_violin_plot
function gamma_violin_plot(sites, posterior_weights, gamma_grid;
    x_label="log(ω + exp(-7))", y_label="Codon Sites",
    sitespace=7.5)

    ypos = [-i * sitespace for i in 1:length(sites)]

    tr(x) = exp(x) - exp(-7) #For consistency, would preferably denote this as trinv
    trinv(x) = log(x + exp(-7))
    tr_prim(x) = exp(x)
    xgrid = gridsetup(exp(-30), exp(7), 100, trinv, tr)
    xinvgrid = trinv.(xgrid)
    point_mass_width = 0.2 #xinvgrid scope

    plot([trinv(1.0), trinv(1.0)], [minimum(ypos) - sitespace, maximum(ypos) + sitespace], color="grey", alpha=0.1, label=:none)

    for i in 1:length(sites)
        posterior_weights_at_site = posterior_weights[i]
        center_line = ypos[i]
        
        # Plot the gamma distributions for the current site
        for j in eachindex(gamma_grid)
            mu, shape, capped = gamma_grid[j]
            weight = posterior_weights_at_site[j]
            G = Gamma(mu, shape)
            y = cappedpdf(G, xgrid, FLAVOR_CAP, point_mass_width, capped=capped, trinv=trinv) .* tr_prim.(xinvgrid) #Account for the change of variables with tr_prim
            #I think we should combine the log+offset transform with a left cap at xgrid[1]
            #This doesn't look quite right: distribute_point_mass!(y, G, xgrid, false, 0.0, point_mass_width, trinv=trinv)
            plot!(xinvgrid, y .+ center_line, fillrange=center_line, alpha=weight, color=:green, linealpha=0, legend=false)
            plot!(xinvgrid, center_line .- y, fillrange=center_line, alpha=weight, color=:green, linealpha=0, legend=false)
        end
    end

    plot!(yticks=(ypos, ["       " for _ in sites]))
    annotate!(0, minimum(ypos) - (sitespace + 2.5), Plots.text(x_label, "black", :center, 10))
    annotate!(minimum(xinvgrid) - 2.5, minimum(ypos) / 2, Plots.text(y_label, "black", :center, 10, rotation=90))

    plot!(ylim=(minimum(ypos) - sitespace, maximum(ypos) + sitespace), xlim = (minimum(xinvgrid), maximum(xinvgrid)))
    
    for i in 1:length(sites)
        annotate!(minimum(xinvgrid) - 0.5, ypos[i], Plots.text("$(sites[i])", "black", :right, 9))
    end
end

#DRY, see collapse_counts
function collapse_weights(param_vec, weight_vec; cases=nothing)
    if isnothing(cases)
        cases = sort(union(param_vec))
    end
    d = Dict(zip(cases, 1:length(cases)))
    storage = zeros(Float64, length(cases))
    for i in 1:length(weight_vec)
        storage[d[param_vec[i]]] += weight_vec[i]
    end
    return storage ./ sum(storage)
end

#Need to make this match the smoothFUBAR setup with the first argument controlling the method (EM, Gibbs, etc)
function FLAVOR(f::FLAVORgrid, outpath; pos_thresh=0.9, verbosity=1, method = (sampler = :DirichletEM, concentration = 0.1, iterations = 2500), plots = true)
    l = size(f.prob_matrix,1)
    num_sites = size(f.prob_matrix,2)
    θ = weightEM(f.prob_matrix, ones(l)./l, conc = method.concentration, iters = method.iterations)

    #Note: we exclude the sites where there were no categories with omega > 1, even if they are uncapped.
    pos_sel_mask = vcat(prop_pos.(f.gridpoints,false),prop_pos.(f.gridpoints,true)) .> 0.0 #[i <= l/2 for i in 1:l]
    pos_prior = sum(pos_sel_mask.*θ)
    posterior_pos = [sum(MolecularEvolution.sum2one(θ .* f.prob_matrix[:,i]).*pos_sel_mask) for i in 1:num_sites];
    bfs = bayes_factor.(posterior_pos,pos_prior)
    detecs_filter = posterior_pos .> pos_thresh
    inv_detecs_filter = posterior_pos .< 1 - pos_thresh

    df = DataFrame()
    df."site" = 1:length(posterior_pos)
    df."P(β>α)" = posterior_pos
    df."BayesFactor" = bfs
    df[!,"P(β>α)>$(pos_thresh)"]= detecs_filter
    df."BayesFactor>10" = bfs .> 10.0
    df."BayesFactor>100" = bfs .> 100.0
    CSV.write(outpath*"_SelectionOutput.csv",df)

    #GammaViolin
    mubyshape = [(g[1], g[2], capped) for capped in [false, true] for g in f.gridpoints]
    cases = [(mu,shape, capped) for capped in [false, true] for mu in f.mugrid for shape in f.shapegrid]
    alloc_grid = θ .* f.prob_matrix

    for (filter, tag) in zip([detecs_filter, inv_detecs_filter], ["pos", "neg"])
        mushape_volumes = Vector{Float64}[]
        sites_to_plot = collect(1:num_sites)[filter]
        for site = sites_to_plot
            mubyshape_weights = collapse_weights(mubyshape, alloc_grid[:, site], cases=cases)
            push!(mushape_volumes, mubyshape_weights)
        end

        if plots
            ysize = 400 + 200 * length(sites_to_plot)
            lmargin = 15 + length(sites_to_plot) / 10
            gamma_violin_plot(sites_to_plot, mushape_volumes, cases)
            plot!(size=(400, ysize), grid=false, left_margin=(lmargin)mm, bottom_margin=10mm)
            savefig(outpath * "_$(tag)_gamma_violin.pdf")
        end
    end

    if plots
        #Posteriors:
        scatter(posterior_pos, label = :none, xlabel = "Codon sites", ylim = (0,1), color=:rainbow_bgyr_35_85_c72_n256, colorbar = :none,
                ylabel = "P(β>α) for some branches", zcolor = posterior_pos .> pos_thresh, markerstrokewidth = 0, size = (800,300), margins = 10Plots.mm)
        for i in 1:length(posterior_pos)
            if posterior_pos[i] > pos_thresh
                annotate!([(i,posterior_pos[i],Plots.text(" $(i)",8,:black,:left))])
            end
        end
        savefig(outpath*"_SitePosteriors.pdf")

        #Bayes factors:
        scatter(bfs, label = :none, xlabel = "Codon sites", color=:rainbow_bgyr_35_85_c72_n256, colorbar = :none,
                ylabel = "Bayes Factor", zcolor = bfs .> 10.0, markerstrokewidth = 0, size = (800,300), margins = 10Plots.mm, yscale = :log10)
        for i in 1:length(bfs)
            if bfs[i] > 10.0
                annotate!([(i,bfs[i],Plots.text(" $(i)",8,:black,:left))])
            end
        end
        savefig(outpath*"_SiteBayesFactors.pdf")
    end

    if verbosity > 0
        for i in 1:length(posterior_pos)
            if posterior_pos[i] > pos_thresh
                println("Site $(i): P(β>α) on some branches = $(round(posterior_pos[i],digits=4))");
            end
        end
    end
    return df
end

export FLAVOR



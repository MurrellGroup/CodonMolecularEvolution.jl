#A call should look like
#f = FAMEgrid(seqnames, seqs, treestring)
#FAME(f, outdir)


struct FAMEgrid{T}
    tr::Function
    trinv::Function
    alphagrid::Vector{T}
    omega1grid::Vector{T}
    omega2grid::Vector{T}
    gridpoints::Vector{Tuple{T,T,T}}
    grid_dims::Vector{Int64}
    prob_matrix::Matrix{T}
    site_scalers::Vector{T}
end

function FAMEgrid(seqnames::Vector{String}, seqs, treestring::String;
    verbosity=1, code=MolecularEvolution.universal_code, optimize_branch_lengths=false)
    tree = FUBAR_init(treestring, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit_2steps(seqnames, seqs, tree, x -> x, code, verbosity=verbosity, optimize_branch_lengths=optimize_branch_lengths)
    
    verbosity > 0 && println("Step 3: Calculating conditional likelihood grid.")

    tr(x) = max(0.0,10^x-0.05)
    trinv(x) =  log10(x+0.05)

    alphagrid = gridsetup(0.01, 10, 8, trinv, tr)
    omega1grid = gridsetup(0.01, 1.0, 14, trinv, tr)
    omega2grid = gridsetup(0.01, 10, 8, trinv, tr)

    #Tiny grid for testing:
    #alphagrid = gridsetup(0.01, 10, 3, trinv, tr)
    #omega1grid = gridsetup(0.01, 1.0, 3, trinv, tr)
    #omega2grid = gridsetup(0.01, 10, 3, trinv, tr)

    #alpha, omega1, omega2
    gridpoints = [(alpha,omega1,omega2) for alpha in alphagrid for omega1 in omega1grid for omega2 in omega2grid];
    grid_dims = length.([alphagrid,omega1grid,omega2grid])

    LL_matrix = marginalized_BWMM_grid(gridpoints,grid_dims, GTRmat, F3x4_freqs, tree, verbosity = verbosity)

    maxi_shift = maximum(LL_matrix,dims = 1)
    prob_matrix = exp.(LL_matrix .- maxi_shift)
    sum_shift = sum(prob_matrix,dims = 1)
    prob_matrix ./= sum_shift
    site_scalers = maxi_shift[:] .+ log.(sum_shift[:])

    return FAMEgrid(tr, trinv, alphagrid, omega1grid, omega2grid, gridpoints, grid_dims, prob_matrix, site_scalers)
end

export FAMEgrid

#Model wrapper - returns models function that gets used in eg. felsenstein_up!
function BWMM_cache(nodelist::Vector{FelNode}, alpha, omega_vec, nuc_matrix, F3x4)
    m = Dict{Float64, BWMModel}()
    for node in nodelist
        haskey(m, node.branchlength) && continue
        Pmatrices = [MolecularEvolution.getPmatrix(GeneralCTMC(MolecularEvolution.MG94_F3x4(alpha,alpha*o,nuc_matrix,F3x4)), node) for o in omega_vec] #We only call getPmatrix once => GeneralCTMC instead of DiagonalizedCTMC
        m[node.branchlength] = BWMModel{PModel}(map(PModel, Pmatrices))
    end
    return m
end

function marginalized_BWMM_grid(gridpoints,grid_dims, GTRmat, F3x4_freqs, tree::FelNode; verbosity = 1, num_parts = 20)
    #Calculates the conditional likelihood grid.
    num_sites = CodonMolecularEvolution.sites(tree.message[1])
    l = prod(grid_dims)
    log_con_lik_matrix = zeros(l,num_sites)

    nodelist = getnodelist(tree)
    branchlengths = unique([node.branchlength for node in nodelist])
    weights = range(0, stop=1, length=num_parts) #Equidistant points between 0 and 1

    for (row_ind, cp) in enumerate(gridpoints)
        alpha, omega1, omega2 = cp
        m = BWMM_cache(nodelist, alpha, [omega1, omega2], GTRmat, F3x4_freqs)
        models = node::FelNode -> [m[node.branchlength]]
        for w in weights #Marginalize over weights
            #Adjust weights in BWMModels
            for branchlength in branchlengths
                m[branchlength].weights = [w, 1-w]
            end
            felsenstein!(tree,models)
            combine!(tree.message[1],tree.parent_message[1])
            log_con_lik_matrix[row_ind,:] .+= MolecularEvolution.site_LLs(tree.message[1]) ./ num_parts
        end
        verbosity > 0 && if mod(row_ind,500)==1
            print(round(100*row_ind/length(gridpoints)),"% ")
            flush(stdout)
        end
    end
    return log_con_lik_matrix
end

function bayes_factor(posterior,prior) #Gonna merge this with FLAVOR regarding DRY
    return (posterior/(1-posterior))*((1-prior)/prior)
end

#Need to make this match the smoothFUBAR setup with the first argument controlling the method (EM, Gibbs, etc)
function FAME(f::FAMEgrid, outpath; pos_thresh=0.9, verbosity=1, method = (sampler = :DirichletEM, concentration = 0.1, iterations = 2500), plots = true)
    l = size(f.prob_matrix,1)
    num_sites = size(f.prob_matrix,2)
    θ = weightEM(f.prob_matrix, ones(l)./l, conc = method.concentration, iters = method.iterations)

    pos_sel_mask = [gridpoint[3] for gridpoint in f.gridpoints] .> 1.0 #Consider letting gridpoints be a matrix (list compr could then be avoided)
    #Downstream is exactly the same as FLAVOR... DRY
    pos_prior = sum(pos_sel_mask.*θ)
    posterior_pos = [sum(MolecularEvolution.sum2one(θ .* f.prob_matrix[:,i]).*pos_sel_mask) for i in 1:num_sites];
    bfs = bayes_factor.(posterior_pos,pos_prior)

    df = DataFrame()
    df."site" = 1:length(posterior_pos)
    df."P(β>α)" = posterior_pos
    df."BayesFactor" = bfs
    df[!,"P(β>α)>$(pos_thresh)"]= posterior_pos .> pos_thresh
    df."BayesFactor>10" = bfs .> 10.0
    df."BayesFactor>100" = bfs .> 100.0
    CSV.write(outpath*"_SelectionOutput.csv",df)

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

export FAME
#TODO: Not limit ourselves to two omega cats
function sim_episodic_div_seqs(alphaD::Distribution, omega1D::Distribution, omega2D::Distribution, weightD::Distribution, numsites::Int, singletree, nucmat::Array{Float64,2}, f3x4::Array{Float64,2};
    scale_total_tree_neutral_expected_subs = -1.0, outpath = "")
    direction = sim_init!(singletree, nucmat, f3x4, scale_total_tree_neutral_expected_subs=scale_total_tree_neutral_expected_subs)

    nucseq_collection = []
    sparams_collection = []
    for i in 1:numsites
        alpha, omegavec, weight = rand(alphaD), [rand(omega1D), rand(omega2D)], rand(weightD)
        betavec = alpha .* omegavec
        weights = Weights([weight, 1-weight])
        models = [DiagonalizedCTMC(MolecularEvolution.MG94_F3x4(alpha, beta, nucmat, f3x4)) for beta in betavec]
        m(n::FelNode) = [sample(models, weights)]
        sample_down!(singletree, m)
        nucs_at_leaves = [n.message[1].obs for n in getleaflist(singletree)]
        push!(nucseq_collection, nucs_at_leaves)
        push!(sparams_collection, (alpha, betavec..., first(weights)))
        lazyprep!(singletree, direction) #This is needed for successive sample_down! calls
    end
    nucseqs = prod(stack(nucseq_collection), dims = 2)[:]
    if outpath != ""
        write_fasta(outpath*".fasta", nucseqs, seq_names = [n.name for n in getleaflist(singletree)])
        open(outpath*".tre", "a") do io write(io, newick(singletree)) end
        write_simparams(outpath, sparams_collection, [:α, :β⁻, :β⁺, :q⁻])
    end
    return nucseqs, [n.name for n in getleaflist(singletree)], singletree, sparams_collection
end

export sim_episodic_div_seqs
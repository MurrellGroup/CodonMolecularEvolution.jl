@doc ("""
    sim_episodic_div_seqs(alphavec::Vector{Float64}, omegamat::Matrix{Float64}, weightmat::Matrix{Float64}, singletree, nucmat::Array{Float64,2}, f3x4::Array{Float64,2};
                            scale_total_tree_neutral_expected_subs = -1.0, outpath = "")

Simulate a set of sequences under a given tree, with a set of alpha values from `alphavec`, omega values from `omegamat`, and probability weights specified in `weightmat`.
For each site `i` and for each branch `b` in the tree, ``P(omega = omegamat[i, j]) = weightmat[i, j]``.
""" * SHARED_SIMDOC)
function sim_episodic_div_seqs(alphavec::Vector{Float64}, omegamat::Matrix{Float64}, weightmat::Matrix{Float64}, singletree, nucmat::Array{Float64,2}, f3x4::Array{Float64,2};
    scale_total_tree_neutral_expected_subs = -1.0, outpath = "")
    @assert size(alphavec, 1) == size(omegamat, 1) == size(weightmat, 1)
    @assert size(omegamat, 2) == size(weightmat, 2)
    ncategories, nnodes = size(omegamat, 2), length(getnodelist(singletree))
    direction = sim_init!(singletree, nucmat, f3x4, scale_total_tree_neutral_expected_subs=scale_total_tree_neutral_expected_subs)

    nucseq_collection = []
    sparams_collection = []
    @views for i in 1:length(alphavec)
        alpha, omegavec, weightvec = alphavec[i], omegamat[i, :], weightmat[i, :]
        models = [DiagonalizedCTMC(MolecularEvolution.MG94_F3x4(alpha, alpha*omega, nucmat, f3x4)) for omega in omegavec]
        model_inds = sample(1:ncategories, Weights(weightvec), nnodes)
        model_ind_dict = Dict(zip(getnodelist(singletree), model_inds))
        m(n::FelNode) = [models[model_ind_dict[n]]]
        sample_down!(singletree, m)
        nucs_at_leaves = [n.message[1].obs for n in getleaflist(singletree)]
        push!(nucseq_collection, nucs_at_leaves)
        push!(sparams_collection, any(i -> omegavec[i] > 1.0 && i ∈ model_inds, 1:ncategories))
        lazyprep!(singletree, direction) #This is needed for successive sample_down! calls
    end
    nucseqs = prod(stack(nucseq_collection), dims = 2)[:]
    if outpath != ""
        write_fasta(outpath*".fasta", nucseqs, seq_names = [n.name for n in getleaflist(singletree)])
        open(outpath*".tre", "a") do io write(io, newick(singletree)) end
        write_simparams(outpath, sparams_collection, [:positive])
    end
    return nucseqs, [n.name for n in getleaflist(singletree)], singletree, sparams_collection
end

export sim_episodic_div_seqs
function sim_episodic_div_seqs(alphaD::Distribution, omegaDs::Vector{<:Distribution}, weightsD::Multinomial, numsites::Int, singletree, nucmat::Array{Float64,2}, f3x4::Array{Float64,2};
    scale_total_tree_neutral_expected_subs = -1.0, outpath = "")
    @assert length(omegaDs) == length(weightsD)
    @assert weightsD.n == length(getnodelist(singletree)) #Modeling-motivated assertion
    direction = sim_init!(singletree, nucmat, f3x4, scale_total_tree_neutral_expected_subs=scale_total_tree_neutral_expected_subs)

    nucseq_collection = []
    sparams_collection = []
    for i in 1:numsites
        alpha, omegavec, weights = rand(alphaD), rand.(omegaDs), rand(weightsD)
        betavec = alpha .* omegavec
        #Fix this multinomial thingy
        models = [DiagonalizedCTMC(MolecularEvolution.MG94_F3x4(alpha, beta, nucmat, f3x4)) for beta in betavec]
        model_ind_dict = Dict(zip(getnodelist(singletree), sample([i for i=1:length(weights) for j=1:weights[i]], weightsD.n, replace=false)))
        m(n::FelNode) = [models[model_ind_dict[n]]]
        sample_down!(singletree, m)
        nucs_at_leaves = [n.message[1].obs for n in getleaflist(singletree)]
        push!(nucseq_collection, nucs_at_leaves)
        push!(sparams_collection, (alpha, betavec..., weights..., any(i -> omegavec[i] > 1.0 && weights[i] > 0, 1:length(omegavec))))
        lazyprep!(singletree, direction) #This is needed for successive sample_down! calls
    end
    nucseqs = prod(stack(nucseq_collection), dims = 2)[:]
    if outpath != ""
        write_fasta(outpath*".fasta", nucseqs, seq_names = [n.name for n in getleaflist(singletree)])
        open(outpath*".tre", "a") do io write(io, newick(singletree)) end
        param_symbols = vcat([:α], Symbol.("β", string.(1:length(omegaDs))), Symbol.("q", string.(1:length(weightsD))), [:positive])
        write_simparams(outpath, sparams_collection, param_symbols)
    end
    return nucseqs, [n.name for n in getleaflist(singletree)], singletree, sparams_collection
end

export sim_episodic_div_seqs
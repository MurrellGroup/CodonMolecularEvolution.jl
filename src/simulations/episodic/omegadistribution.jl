function sim_episodic_div_seqs(alphaD::Distribution, omegaD_template::T, param1D::UnivariateDistribution, param2D::UnivariateDistribution, numsites::Int, singletree, nucmat::Array{Float64,2}, f3x4::Array{Float64,2};
    scale_total_tree_neutral_expected_subs = -1.0, outpath = "") where {T<:UnivariateDistribution}
    @assert length(params(omegaD_template)) == 2 "The omegaD must have two degrees of freedom (two parameters), one accounting for α and the other β"
    direction = sim_init!(singletree, nucmat, f3x4, scale_total_tree_neutral_expected_subs=scale_total_tree_neutral_expected_subs)

    nucseq_collection = []
    sparams_collection = []
    for i in 1:numsites
        alpha, param1, param2 = rand(alphaD), rand(param1D), rand(param2D)
        omegaD = T(param1, param2)
        omegas = rand(omegaD, length(getnodelist(singletree)))
        omega_dict = Dict(zip(getnodelist(singletree), omegas))
        m(n::FelNode) = [GeneralCTMC(MolecularEvolution.MG94_F3x4(alpha, alpha*omega_dict[n], nucmat, f3x4))]

        #Need some test to determine if the site is a positive or not
        #Draw omegas prior to the sample_down! call and cache them, and use the cache for m func

        sample_down!(singletree, m)
        nucs_at_leaves = [n.message[1].obs for n in getleaflist(singletree)]
        push!(nucseq_collection, nucs_at_leaves)
        push!(sparams_collection, (alpha, (param1, param2), any(omega -> omega > 1.0, omegas)))
        lazyprep!(singletree, direction) #This is needed for successive sample_down! calls
    end
    nucseqs = prod(stack(nucseq_collection), dims = 2)[:]
    if outpath != ""
        write_fasta(outpath*".fasta", nucseqs, seq_names = [n.name for n in getleaflist(singletree)])
        open(outpath*".tre", "a") do io write(io, newick(singletree)) end
        write_simparams(outpath, sparams_collection, [:α, Symbol(T), :positive])
    end
    return nucseqs, [n.name for n in getleaflist(singletree)], singletree, sparams_collection
end

export sim_episodic_div_seqs
#TODO: DRY, see sim_alphabeta_seqs
function sim_episodic_div_seqs(alphaD::Distribution, muD::Distribution, shapeD::Distribution, numsites::Int, singletree, nucmat::Array{Float64,2}, f3x4::Array{Float64,2};
    scale_total_tree_neutral_expected_subs = -1.0, outpath = "")
    eq_freqs = MolecularEvolution.F3x4_eq_freqs(f3x4)
    eq_partition = CodonPartition(1)
    eq_partition.state .= eq_freqs
    lazy_template = LazyPartition{CodonPartition}()
    internal_message_init!(singletree, lazy_template)
    direction = LazyDown(isleafnode)
    lazyprep!(singletree, eq_partition, direction=direction)

    if scale_total_tree_neutral_expected_subs > 0.0
        total_bl = sum([n.branchlength for n in getnodelist(singletree)])
        testmat = MolecularEvolution.MG94_F3x4(1.0, 1.0, nucmat, f3x4)
        current_expected_neutral_subs = (-sum(eq_freqs .* diag(testmat)))*total_bl
        scale = scale_total_tree_neutral_expected_subs/current_expected_neutral_subs
        for n in getnodelist(singletree)
            n.branchlength *= scale
        end
    end

    nucseq_collection = []
    cparams_collection = []
    for i in 1:numsites
        alpha, mu, shape = rand(alphaD), rand(muD), rand(shapeD)
        omegaD = Gamma(mu, shape)
        omegas = rand(omegaD, length(getnodelist(singletree)))
        omega_dict = Dict(zip(getnodelist(singletree), omegas))
        m(n::FelNode) = [GeneralCTMC(MolecularEvolution.MG94_F3x4(alpha, alpha*omega_dict[n], nucmat, f3x4))]

        #Need some architecture to determine if the site is a positive or not
        #Draw omegas prior to the sample_down! call and cache them, and use the cache for m func

        sample_down!(singletree, m)
        nucs_at_leaves = [n.message[1].obs for n in getleaflist(singletree)]
        push!(nucseq_collection, nucs_at_leaves)
        push!(cparams_collection, (alpha, mu, shape, any(omega -> omega > 1.0, omegas)))
        lazyprep!(singletree, direction) #This is needed for successive sample_down! calls
    end
    nucseqs = prod(stack(nucseq_collection), dims = 2)[:]
    if outpath != ""
        write_fasta(outpath*".fasta", nucseqs, seq_names = [n.name for n in getleaflist(singletree)])
        open(outpath*".tre", "a") do io write(io, newick(singletree)) end
        df = DataFrame()
        alphas = [cp[1] for cp in cparams_collection]
        df."site" = 1:numsites
        df."α" = alphas
        df."mu" = [cp[2] for cp in cparams_collection]
        df."shape" = [cp[3] for cp in cparams_collection]
        df."positive" = [cp[4] for cp in cparams_collection]
        CSV.write(outpath*"_SimParams.csv",df)
    end
    return nucseqs, [n.name for n in getleaflist(singletree)], singletree, cparams_collection
end

export sim_episodic_div_seqs
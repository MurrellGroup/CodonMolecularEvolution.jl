#Estimates a tree and nuc model from data, calibrated so that alpha=1
function alphabeta_setup(seqnames, seqs, treestring;
    verbosity=1, code=MolecularEvolution.universal_code, optimize_branch_lengths=false)
    tree, _ = FUBAR_init("", treestring, exports=false, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit_2steps(seqnames, seqs, tree, x -> x, code, verbosity=verbosity, optimize_branch_lengths=optimize_branch_lengths)
    return tree, GTRmat, F3x4_freqs
end

export alphabeta_setup

"""
    sim_alphabeta_seqs(alphavec::Vector{Float64}, betavec::Vector{Float64}, singletree, nucmat::Array{Float64,2}, f3x4::Array{Float64,2};
                            scale_total_tree_neutral_expected_subs = -1.0, outpath = "")

Simulate a set of sequences under a given tree, with a set of alpha and beta values.
f3x4 is a 3-by-4 matrix of position-specific nucleotide frequencies.
nucmat is a 4-by-4 matrix of nucleotide substitution rates.
If scale_total_tree_neutral_expected_subs > 0, then the tree is scaled so that if alpha=beta=1 for all sites, the expected number of neutral substitutions is equal to scale_total_tree_neutral_expected_subs.
The sequences are written to a fasta file, and the tree is written to a newick file.
"""
function sim_alphabeta_seqs(alphavec::Vector{Float64}, betavec::Vector{Float64}, singletree, nucmat::Array{Float64,2}, f3x4::Array{Float64,2};
                            scale_total_tree_neutral_expected_subs = -1.0, outpath = "")
    @assert length(alphavec) == length(betavec)
    direction = sim_init!(singletree, nucmat, f3x4, scale_total_tree_neutral_expected_subs=scale_total_tree_neutral_expected_subs)

    nucseq_collection = []
    for i in 1:length(alphavec)
        alpha,beta = alphavec[i], betavec[i]
        m = DiagonalizedCTMC(MolecularEvolution.MG94_F3x4(alpha, beta, nucmat, f3x4))
        sample_down!(singletree, m)
        nucs_at_leaves = [n.message[1].obs for n in getleaflist(singletree)]
        push!(nucseq_collection, nucs_at_leaves)
        lazyprep!(singletree, direction) #This is needed for successive sample_down! calls
    end
    nucseqs = prod(stack(nucseq_collection), dims = 2)[:]
    if outpath != ""
        write_fasta(outpath*".fasta", nucseqs, seq_names = [n.name for n in getleaflist(singletree)])
        open(outpath*".tre", "a") do io write(io, newick(singletree)) end
    end
    return nucseqs, [n.name for n in getleaflist(singletree)], singletree
end

export sim_alphabeta_seqs

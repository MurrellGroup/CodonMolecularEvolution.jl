#Estimates a tree and nuc model from data, calibrated so that alpha=1
function alphabeta_setup(seqnames, seqs, treestring;
    verbosity=1, code=MolecularEvolution.universal_code, optimize_branch_lengths=false)
    tree, _ = FUBAR_init("", treestring, exports=false, verbosity=verbosity)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit_2steps(seqnames, seqs, tree, x -> x, code, verbosity=verbosity, optimize_branch_lengths=optimize_branch_lengths)
    return tree, GTRmat, F3x4_freqs
end

export alphabeta_setup

#These were derived from a Flu dataset
const demo_f3x4 = [ 0.293117 0.184379 0.295274 0.190878;
                    0.342317 0.199907 0.154328 0.267101;
                    0.231987 0.217801 0.241637 0.272234]


const demo_nucmat = [  -0.256236 0.0697056 0.152411 0.034119;
                        0.0697056 -0.274119 0.0596187 0.144795;
                        0.152411 0.0596187 -0.251381 0.0393506;
                        0.034119 0.144795 0.0393506 -0.218264]


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

include("alphabeta/alphabeta.jl")
include("episodic/omegamixture.jl")
include("episodic/omegadistribution.jl")

#These were derived from a Flu dataset
const demo_f3x4 = [ 0.293117 0.184379 0.295274 0.190878;
                    0.342317 0.199907 0.154328 0.267101;
                    0.231987 0.217801 0.241637 0.272234]


const demo_nucmat = [  -0.256236 0.0697056 0.152411 0.034119;
                        0.0697056 -0.274119 0.0596187 0.144795;
                        0.152411 0.0596187 -0.251381 0.0393506;
                        0.034119 0.144795 0.0393506 -0.218264]

#Should maybe be in MolecularEvolution.jl
function standard_tree_sim(ntaxa)
    n(t) = (10*ntaxa)/(1+exp(t-10))
    return sim_tree(ntaxa,n,ntaxa/5, mutation_rate = 0.05)
end
function ladder_tree_sim(ntaxa)
    n(t) = ntaxa/10
    return sim_tree(ntaxa,n,1.0, mutation_rate = 1/sqrt(ntaxa))
end

export standard_tree_sim, ladder_tree_sim

function sim_init!(singletree, nucmat::Array{Float64}, f3x4::Array{Float64}; scale_total_tree_neutral_expected_subs = -1.0)
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
        rescale_branchlengths!(singletree, scale)
    end
    return direction
end

function write_simparams(outpath, sparams_collection, colnames::Vector{Symbol})
    df = DataFrame()
    numsites = length(sparams_collection)
    df."site" = 1:numsites

    for (i, colname) in enumerate(colnames)
        df[!, colname] = [sp[i] for sp in sparams_collection]
    end

    CSV.write(outpath * "_SimParams.csv", df)
end
include("alphabeta/alphabeta.jl")
include("episodic/omegamixture.jl")
include("episodic/omegadistribution.jl")

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
        for n in getnodelist(singletree)
            n.branchlength *= scale
        end
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
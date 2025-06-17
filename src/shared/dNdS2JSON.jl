# We assume order doesn't matter so we just use `Dict`, otherwise `DataStructures.OrderedDict` could be used.
# I assume there's only one partition (HyPhy) "0" (can be generalized later)

const SHAREDJSONFIELDS = Dict(
    "analysis" => ("info", "version", "citation", "authors", "contact", "requirements"),
    "input" => ("file name", "number of sequences", "number of sites", "partition count", "trees"),
    "fits" => ("Log Likelihood", "estimated parameters", "AIC-c", "Rate Distributions", "Nucleotide GTR"),
    "data partitions" => ("0",),
    "branch attributes" => ("0",),
    "tested" => nothing,  # specified later
    "timers" => nothing  # specified later
)

#=
Must contain the field `d`
A NamedTuple with (atleast) the fields: seqs, seqnames, treestring, timers should be passed along downstream
Must contain the method `model_stages(::Type{ConcretedNdS2JSON})` (should be static)
=#
abstract type dNdS2JSON end

# Rely on input2JSON! already having been called
sites(method::dNdS2JSON) = method.d["input"]["number of sites"]

#= 
FUBAR-esque dNdS methods
Let us denote a concrete type by `ConcretedNdS2JSON <: FUBAREsque2JSON`
Must contain the method `expanded_headers(::Type{ConcretedNdS2JSON})`
A NamedTuple with (atleast) the fields: df, θ, posterior_mat, categories should be passed along downstream
=#
abstract type FUBAREsque2JSON <: dNdS2JSON end

struct difFUBAR2JSON <: FUBAREsque2JSON
    d::Dict{String, Any}
    model_stages::Tuple{Vararg{String}}

    function difFUBAR2JSON()
        new(Dict{String, Any}())
        #We might need this... checking with Steve
        #new(Dict{String, Any}("analysis" => Dict{String, Any}(zip(SHAREDJSONFIELDS["analysis"], analysis_info(difFUBAR2JSON)))))
    end
end

analysis_info(::Type{difFUBAR2JSON}) = (
    "Perform a site-wise comparison of evolutionary pressure between two selected sets of branches.",
    "version...",
    "citation...",
    "authors...",
    "contact...",
    "requirements..."
)

expanded_headers(::Type{difFUBAR2JSON}) = (
    "The codon site of interest",
    "Posterior probability of dNdS for branch group 1 > dNdS for branch group 2",
    "Posterior probability of dNdS for branch group 2 > dNdS for branch group 1",
    "Posterior probability of positive selection for branch group 1 at a site",
    "Posterior probability of positive selection for branch group 2 at a site",
    "Mean posterior synonymous substitution rate at a site",
    "Mean posterior dNdS for branch group 1 at a site",
    "Mean posterior dNdS for branch group 2 at a site"
)

model_stages(::Type{difFUBAR2JSON}) = (
    "Overall",
    "Global Fit",
    "Grid Calculations",
    "MCMC on grid weight allocations"
)

function dNdS2JSON!(method::dNdS2JSON, nt::NamedTuple)
    shared2JSON!(method, nt)
    specialized2JSON!(method, nt)
end

function dNdS2JSON(method::dNdS2JSON, nt::NamedTuple)
    dNdS2JSON!(method, nt)
    json = JSON.json(method.d, 2)
    if nt.exports
        open(nt.outpath * ".json", "w") do io
            write(io, json)
        end
    end
    return json
end

function input2JSON!(method::dNdS2JSON, nt::NamedTuple)
    seqnames, seqs, treestring = nt.seqnames, nt.seqs, nt.treestring
    input_info = ("<CodonMolEv difFUBAR doesn't know>", length(seqnames), Int(length(seqs[1]) / 3), 1, Dict("0" => treestring))
    method.d["input"] = Dict{String, Any}(zip(SHAREDJSONFIELDS["input"], input_info))
end

# 1 partition as default/fallback
function datapartitions2JSON!(method::dNdS2JSON, nt::NamedTuple)
    num_sites = sites(method)
    method.d["data partitions"] = Dict{String, Any}(
        "0" => Dict{String, Any}(
            "name" => "default",
            "sites" => collect(1:num_sites)
            )
        )
end

function tested2JSON!(method::dNdS2JSON, nt::NamedTuple)
    method.d["tested"] = 0
end

function timers2JSON!(method::dNdS2JSON, nt::NamedTuple)
    inner_dict = Dict{String, Any}()
    for (i, (stage, timer)) in enumerate(zip(model_stages(typeof(method)), nt.timers))
        inner_dict[stage] = Dict("timer" => timer, "order" => i-1)
    end
    method.d["timers"] = inner_dict
end

function shared2JSON!(method::dNdS2JSON, nt::NamedTuple)
    input2JSON!(method, nt)
    fits2JSON!(method, nt)
    datapartitions2JSON!(method, nt)
    branchattributes2JSON!(method, nt)
    tested2JSON!(method, nt)
    timers2JSON!(method, nt)
end

function canonical_specialized2JSON!(method::FUBAREsque2JSON, nt::NamedTuple)
    df, θ, posterior_mat, categories = nt.df, nt.θ, nt.posterior_mat, nt.categories
    # MLE
    headers = Tuple(zip(names(df), expanded_headers(typeof(method))))
    method.d["MLE"] = Dict("headers" => headers, "content" => Dict("0" => eachrow(Matrix(df))))

    # Grid
    method.d["grid"] = eachrow(hcat(categories, θ))

    # Posterior
    method.d["posterior"] = Dict("0" => Dict(zip(0:sites(method)-1, eachcol(posterior_mat))))
end

function specialized2JSON!(method::difFUBAR2JSON, nt::NamedTuple)
    canonical_specialized2JSON!(method, nt)
    # Settings
    method.d["settings"] = Dict{String, Any}(
        "chain-length" => nt.iters,
        "burn-in" => nt.burnin,
        "concentration" => nt.concentration,
        "posterior" => nt.pos_thresh
    )
end

function fits2JSON!(method::difFUBAR2JSON, nt::NamedTuple)
    method.d["fits"] = Dict{String, Any}("Log Likelihood" => nt.LL)
end


function branchattributes2JSON!(method::difFUBAR2JSON, nt::NamedTuple)
    # Branch attributes
    node_dict = Dict{String, Any}()
    for n in nodes(nt.tree)
        nd = Dict{String, Any}("Nucleotide GTR" => n.branchlength)
        model_ind = CodonMolecularEvolution.model_ind(n.name, nt.tags)
        nd["Branch group"] = model_ind > length(nt.tags) ? "background" : string(model_ind)
        if isleafnode(n)
            nd["original name"] = nt.leaf_name_transform(n.name)
        end
        node_dict[n.name] = nd
    end
    # Attributes description
    attr_dict = Dict{String, Any}(
        "Nucleotide GTR" => Dict("attribute type" => "branch length", "display order" => 0),
        "Branch group" => Dict("attribute type" => "node label", "display order" => 1),
        "original name" => Dict("attribute type" => "node label", "display order" => -1)
    )
    method.d["branch attributes"] = Dict("0" => node_dict, "attributes" => attr_dict)
end

function tested2JSON!(method::difFUBAR2JSON, nt::NamedTuple)
    tags = nt.tags
    f(s::String) = ifelse(CodonMolecularEvolution.model_ind(s, tags) > length(tags), "background", "test")
    n = map(x -> x.name, nodes(nt.tree))
    method.d["tested"] = Dict("0" => Dict(zip(n, f.(n))))
end

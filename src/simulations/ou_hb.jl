#Halpern and Bruno model where amino acid fitnesses evolve over time using a piecewise constant approximation to an OU process.
#Authors: Hassan Sadiq and Ben Murrell

#To do:
#- Check that an alternate genetic code will thread through all of this.
#- Introduce a jump model with an additional class of jumps to independent draws from the equilibrium fits, and make it easy to use this instead
#- Allow viz functions to have both kids of offsets (include the AA offsets in the fitness plot, but not the codon offsets)

#TODO fix t < time + next_event double counting issue

#TODO fix pushing into coll within if statement

abstract type CodonSimulationPartition <: Partition end #Must have .sites::Int64, .codons::Vector{Int64}, and .code::GeneticCode

#Does nothing because the `forward!` function implicitly samples.
function MolecularEvolution.sample_partition!(p::CodonSimulationPartition)
end

function MolecularEvolution.partition2obs(part::CodonSimulationPartition)
    return join([part.code.sense_codons[part.codons[i]] for i in 1:part.sites])
end

const d_code = MolecularEvolution.universal_code

"""
    HB_fixation_rate(from_codon, to_codon)

    Returns the fixation rate of a mutation from `from_codon` to `to_codon` under the HB98 model.
"""
function HB_fixation_rate(from_codon, to_codon)
    diff = to_codon - from_codon
    if abs(diff) < 0.001
        f_ab = 1/12 * (12 + diff*(6 + diff))
    else
        f_ab = diff/(1-exp(-diff))
    end
end

#=
#Confirming local approx
t0(diff) = diff/(1-exp(-diff))
t1(diff) = 1/12 * (12 + diff*(6 + diff))
r = -0.01:0.000001:0.01
plot(r, t0.(r) .- t1.(r), label = "1/12 * (12 + diff*(6 + diff))")
=#

"""
    HB98AA_matrix(alpha, nuc_matrix, AA_fitness; genetic_code = MolecularEvolution.universal_code)

Returns the rate matrix for a codon model using the HB98 model where each AA has a different fitness. `alpha` is the synonymous rate, `nuc_matrix` is the symmetric nucleotide substitution matrix,
and `AA_fitness` is the fitness of each amino acid.
"""
function HB98AA_matrix(alpha, nuc_matrix, AA_fitness; genetic_code = MolecularEvolution.universal_code)
    codon_matrix = zeros(length(genetic_code.sense_codons), length(genetic_code.sense_codons))
    for p in genetic_code.syn_positions
        codon_matrix[p[1][1], p[1][2]] = alpha * nuc_matrix[p[2][2], p[2][3]]
    end
    for p in genetic_code.nonsyn_positions
        f_ab = HB_fixation_rate(AA_fitness[genetic_code.codon2AA_pos[p[1][1]]], AA_fitness[genetic_code.codon2AA_pos[p[1][2]]])
        codon_matrix[p[1][1], p[1][2]] = alpha * nuc_matrix[p[2][2], p[2][3]] * f_ab
    end
    for i = 1:length(genetic_code.sense_codons)
        codon_matrix[i, i] = -sum(codon_matrix[i, :])
    end
    return codon_matrix
end

#Assumption: nuc matrix is symmetric
function HB98_eqfreqs(fitnesses) #Result from HB98
    v = exp.(fitnesses .- maximum(fitnesses))
    return v ./ sum(v)
end
HB98AA_eqfreqs(AA_fitness; genetic_code = MolecularEvolution.universal_code) = HB98_eqfreqs(AA_fitness[genetic_code.codon2AA_pos])
function HB98AA_expected_subs_per_site(nuc_matrix, AA_fitness; genetic_code = MolecularEvolution.universal_code)
    eqfreqs = HB98AA_eqfreqs(AA_fitness, genetic_code = genetic_code)
    Q = HB98AA_matrix(1.0, nuc_matrix, AA_fitness, genetic_code = genetic_code)
    return expected_subs_per_site(Q, eqfreqs)
end
HB98AA_neutral_scale_nuc_matrix(nucm; genetic_code = MolecularEvolution.universal_code) = nucm ./= HB98AA_expected_subs_per_site(nucm, zeros(20), genetic_code = genetic_code)
#/Assumption

"""
    PiecewiseOUModel(event_rate::Float64, eq_std::Float64, mixing::Float64; delta_t = 1.0)
    PiecewiseOUModel(offsets::Vector{Float64})
    PiecewiseOUModel(event_rate::Float64, eq_std::Float64, mixing::Float64, mu::Union{Float64,Vector{Float64}}, offsets::Union{Float64,Vector{Float64}}, codon_offsets::Union{Float64,Vector{Float64}})

A piecewise constant approximation to an OU process, intended to simulate fitnesses evolving over phylogenies.
The equilibrium standard deviation is directly parameterized (`eq_std`), as is the rate at which the process mixes to equilibrium (`mixing`).
`event_rate` controls how often the fitness changes occur, where the mixing rate is scaled to compensate for the increased rate of change to achieve
approximately the same amount of change per unit time even as the `event_rate` changes. A very high `event_rate` will behave more like continuous diffusion,
but will be more computationally expensive to sample from. `mu` can also be set to control the mean fitnesses.
The model also permits `offsets`, which are added to the fitnesses as they are passed into the model. For a single process, these are confounded with the mean `mu`
but if the offsets change (eg. from one branch to another) the effective fitnesses will immidiately change, whereas if `mu` changes the fitnesses will drift towards `mu`.
"""
mutable struct PiecewiseOUModel{M<:Union{Float64,AbstractVector{Float64}},O<:Union{Float64,AbstractVector{Float64}},C<:Union{Float64,AbstractVector{Float64}}} 
    mu::M #OU mean. Scalar (in which case it is the same for all AAs) or vector (one per AA)
    delta_t::Float64 #Scales the size of the fitness jumps
    var::Float64 #OU instantaneous diffusion variance
    reversion_rate::Float64 #Rate of reversion to the mean
    event_rate::Float64 #Rate of fitness jumps
    offsets::O #These do not change, thus allowing preferred codons. Incorporated right at the end.
    eq_std::Float64 #Standard deviation of the equilibrium distribution
    codon_offsets::C #These do not change, thus allowing preferred codons.
    function PiecewiseOUModel(event_rate::Float64, eq_std::Float64, mixing::Float64; delta_t = 1.0)
        new{Float64,Float64,Float64}(0.0, delta_t, 2*(mixing/event_rate)*(eq_std^2), mixing/event_rate, event_rate, 0.0, eq_std, 0.0)
    end
    #Need to test this constructor:
    function PiecewiseOUModel(event_rate::Float64, eq_std::Float64, mixing::Float64, mu::M, offsets::O, codon_offsets::C) where {M,O,C}
        new{M,O,C}(mu, 1.0, 2*(mixing/event_rate)*(eq_std^2), mixing/event_rate, event_rate, offsets, eq_std, codon_offsets)
    end
    function PiecewiseOUModel(offsets::Vector{Float64}) #Constructor for static model, using just offsets for fitnesses
        new{Float64,Vector{Float64}, Float64}(0.0, 1.0, 1e-15, 1.0, 0.0, offsets, 1e-15, 0.0)
    end
end

"""
    jump(x_source, m::PiecewiseOUModel)

Evolves values over time using a piecewise constant approximation to an OU process, where this function computes the new distribution for a single discrete jump.
`x_source` is the vector of fitnesses, and m is the PiecewiseOUModel.

"""
function jump(x_source, x_rand, m::PiecewiseOUModel)
    a = (m.reversion_rate * m.delta_t)
    return x_rand .* (m.eq_std * sqrt(1 - exp(-2*a))) .+ m.mu .+ ((x_source .- m.mu) .* exp(-a))
end
jump(x_source::Float64, m::PiecewiseOUModel) = jump(x_source, randn(), m)
jump(x_source::Vector{Float64}, m::PiecewiseOUModel) = jump(x_source, randn(length(x_source)), m)

#=
m = PiecewiseOUModel(100.0, 10.0, 10.0)
CodonMolecularEvolution.jump(zeros(20), m)
=#

"""
    HB98AA_row(current_codon, alpha, nuc_matrix, AA_fitness; genetic_code=MolecularEvolution.universal_code)

Returns the rate row for a codon model using the HB98 model where each AA has a different fitness. `current_codon` is the current codon, `alpha` is the synonymous rate,
`nuc_matrix` is the symmetric nucleotide substitution matrix, and `AA_fitness` is the fitness of each amino acid.
"""
function HB98_row(current_codon, alpha, nuc_matrix, fitness; genetic_code=MolecularEvolution.universal_code)
    row = zeros(genetic_code.num_sense)
    codon_aa_i = fitness[current_codon]
    for j in 1:genetic_code.num_sense
        if genetic_code.codon_to_nuc_map[current_codon,j] != (-1,-1,-1)
            c2n_map = genetic_code.codon_to_nuc_map[current_codon,j]
            f_ab = HB_fixation_rate(codon_aa_i,fitness[j])
            row[j] = alpha * nuc_matrix[c2n_map[2],c2n_map[3]] * f_ab
        end
    end
    return row
end

HB98AA_row(current_codon, alpha, nuc_matrix, AA_fitness, codon_fitness_offsets; genetic_code=MolecularEvolution.universal_code) = HB98_row(current_codon, alpha, nuc_matrix, AA_fitness[genetic_code.codon2AA_pos] .+ codon_fitness_offsets, genetic_code = genetic_code)

#= Check the row matches the matrix:
fs = randn(20)
c = rand(1:61)
r1 = CodonMolecularEvolution.HB98AA_row(c, 1.23, CodonMolecularEvolution.demo_nucmat, fs, 0.0)
q1 = CodonMolecularEvolution.HB98AA_matrix(1.23, CodonMolecularEvolution.demo_nucmat, fs)
q1[c,c] = 0
@assert isapprox(sum(q1[c,:] .- r1), 0.0)
=#

"""
    jumpy_HB_codon_evolve(fitnesses, codon, ou_model, nuc_matrix, alpha, time;
        genetic_code = MolecularEvolution.universal_code, push_into = nothing)

Evolves fitnesses and codons over time using the HB98 model. `fitnesses` is the vector of fitnesses, `codon` is the current codon, `ou_model` is the OU model,
`nuc_matrix` is the symmetric nucleotide substitution matrix, `alpha` is the synonymous rate, and `time` is the total time to evolve over.
"""
function jumpy_HB_codon_evolve(fitnesses, codon, scaled_fitness_model, nuc_matrix, alpha, time;
        genetic_code = MolecularEvolution.universal_code, push_into = nothing, logNe = nothing, time_origin = 0.0)
    if isnothing(logNe)
        transform = identity
    else
        transform = x -> x .* exp(logNe) .* 2
    end
    codon_jumps = 0
    fitness_jumps = 0
    current_fits = copy(fitnesses)
    current_codon = codon
    t = 0.0
    next_event = 0.0
    while t+next_event < time
        HBrow = HB98AA_row(current_codon, alpha, nuc_matrix, transform(current_fits .+ scaled_fitness_model.offsets), transform(scaled_fitness_model.codon_offsets), genetic_code=genetic_code)
        sum_HBrow = sum(HBrow)
        rOU,rHB = (scaled_fitness_model.event_rate,sum_HBrow)
        total_rate = rOU+rHB
        next_event = randexp()/total_rate
        t = t+next_event
        if t < time
            event_index = sample(1:2,Weights([rOU,rHB])) 
            if event_index == 1 # Fitness jump event
                fitness_jumps += 1
                current_fits = jump(current_fits, scaled_fitness_model)
            else # Codon substitution event
                codon_jumps += 1
                current_codon = sample(1:length(HBrow),Weights(HBrow))
            end
        end
        if !isnothing(push_into)
            if !isnothing(logNe)
                push!(push_into,(time_origin+t,current_codon,transform(copy(current_fits)),logNe))
            else
                push!(push_into,(time_origin+t,current_codon,transform(copy(current_fits))))
            end
        end
    end
    return current_fits, current_codon, codon_jumps, fitness_jumps
end


"""
    ShiftingHBSimModel(sites, alphas, ou_params, nuc_matrix; rescale = true)

A model for simulating fitnesses evolving over phylogenies using the HB98 model. `sites` is the number of sites, `alphas` is a vector of synonymous rates (one per site),
`ou_params` is a vector of `PiecewiseOUModel`s (one per site), and `nuc_matrix` is the symmetric nucleotide substitution matrix (shared across sites).
If 'rescale' is true, then the nuc matrix is scaled so that, when `alpha=1` and the fitnesses`=0`, the HB model expects one substitution per site per unit time.
"""
mutable struct ShiftingHBSimModel <: MolecularEvolution.SimulationModel
    sites::Int64
    alphas::Vector{Float64}
    ou_params::Vector{PiecewiseOUModel}
    nuc_matrix::Matrix{Float64}
    code::MolecularEvolution.GeneticCode
    function ShiftingHBSimModel(s,a,oup,n; rescale = true, code = MolecularEvolution.universal_code)
        new(s,a,oup,rescale ? HB98AA_neutral_scale_nuc_matrix(n, genetic_code = code) : n, code)
    end
end

ShiftingHBSimModel(nuc_matrix::Matrix{Float64}, ou_params::Vector{PiecewiseOUModel{A,B,C}}; alpha = ones(length(ou_params)), rescale = true, code = MolecularEvolution.universal_code) where {A,B,C} = ShiftingHBSimModel(length(ou_params), alpha, ou_params, nuc_matrix, rescale = rescale, code = code)


"""
    ShiftingHBSimPartition(model::ShiftingHBSimModel; burnin_time = 100.0, code = MolecularEvolution.universal_code)

Constructs a partition that tracks evolving fitnesses and codons. Only useable for sampling (not likelihood calculations).
"""
mutable struct ShiftingHBSimPartition <: CodonSimulationPartition
    sites::Int64
    fitnesses::Matrix{Float64}
    codons::Vector{Int64}
    code::MolecularEvolution.GeneticCode
    function ShiftingHBSimPartition(sites; code = MolecularEvolution.universal_code)
        new(sites,zeros(length(code.amino_acids),sites),ones(Int64,sites), code)
    end
    function ShiftingHBSimPartition(model::ShiftingHBSimModel; burnin_time = 100.0)
        fits = zeros(length(model.code.amino_acids),length(model.ou_params))
        for (i,m) in enumerate(model.ou_params)
            fits[:,i] .= randn(length(model.code.amino_acids)) .* m.eq_std .+ m.mu
        end
        #Starting codons ignore codon_offsets. Should be ok because burn-in, but should adjust.
        codons = [sample(1:61, Weights(HB98AA_eqfreqs(fits[:,i] .+ model.ou_params[i].codon_offsets))) for i in 1:length(model.ou_params)]
        for (i,m) in enumerate(model.ou_params)
            f, c, _, _ = jumpy_HB_codon_evolve(fits[:,i], codons[i], m, model.nuc_matrix, model.alphas[i], burnin_time, genetic_code = model.code)
            fits[:,i] .= f
            codons[i] = c
        end
        new(length(model.ou_params), fits, codons, model.code)
    end
end


function MolecularEvolution.forward!(dest::ShiftingHBSimPartition,
        source::ShiftingHBSimPartition,
        model::ShiftingHBSimModel,
        node::FelNode)
    for site in 1:model.sites
        fitnesses = source.fitnesses[:,site]
        codon = source.codons[site]
        ou_model = model.ou_params[site]
        alpha = model.alphas[site]
        fitnesses,codon,_,_ = jumpy_HB_codon_evolve(fitnesses,
            codon,
            ou_model,
            model.nuc_matrix,
            alpha,
            node.branchlength,
            genetic_code = model.code)
        dest.fitnesses[:,site] .= fitnesses
        dest.codons[site] = codon
    end
end

"""
    dNdS(q1, q0, p; code = MolecularEvolution.universal_code)
    dNdS(q1, q0; code = MolecularEvolution.universal_code)

Returns an analytic expectation of the dN/dS ratio for a codon model. `q1` is the rate matrix where selection is active (eg. a Halpern and Bruno model with a set of fitnesses),
and `q0` is the corresponding rate matrix where selection is inactive (eg. a Halpern and Bruno model with all fitnesses equal).
`p` is the frequency distribution over codons that the dN/dS ratio is computed against. If not provided, this is computed as the equilibrium from `q1`.
If only a vector of fitnesses are provided, then the `q1`, `q0`, and `p` are computed assuming a Halpern and Bruno model.
```
fs = randn(20)
nucm = CodonMolecularEvolution.demo_nucmat
q1 = HB98AA_matrix(1.0, nucm, fs)
q0 = HB98AA_matrix(1.0, nucm, zeros(20))
dNdS(q1, q0)
```
"""
function dNdS(q1, q0, p::Vector{Float64}; code = d_code)
    nsi = [CartesianIndex(i[1][1], i[1][2]) for i in code.nonsyn_positions]
    si = [CartesianIndex(i[1][1], i[1][2]) for i in code.syn_positions]
    dN_numerator = sum((q1 .* p)[nsi])
    dN_denominator = sum((q0 .* p)[nsi])
    dS_numerator = sum((q1 .* p)[si])
    dS_denominator = sum((q0 .* p)[si])
    dn = (dN_numerator/dN_denominator)
    ds = (dS_numerator/dS_denominator)
    return dn/ds
end

dNdS(q1, q0; code = MolecularEvolution.universal_code) = dNdS(q1, q0, exp(q1 * 100)[1,:], code = code)

"""
    HBdNdS(fs::Vector{Float64}; code = MolecularEvolution.universal_code, nucm = CodonMolecularEvolution.demo_nucmat)
    HBdNdS(fs_pre::Vector{Float64}, fs_post::Vector{Float64}; code = MolecularEvolution.universal_code, nucm = CodonMolecularEvolution.demo_nucmat)

Returns the expected dN/dS ratio for a Halpern and Bruno model with a vector of fitnesses. If two vectors are provided, then the dN/dS ratio is computed for the shift from `fs_pre` to `fs_post`.
"""
function HBdNdS(fs::Vector{Float64}; code = d_code, nucm = CodonMolecularEvolution.demo_nucmat) 
    q1 = HB98AA_matrix(1.0, nucm, fs, genetic_code = code)
    q0 = HB98AA_matrix(1.0, nucm, zeros(20), genetic_code = code)
    return dNdS(q1, q0, HB98AA_eqfreqs(fs, genetic_code = code), code = code)
end

HBdNdS(fs_pre::Vector{Float64}, fs_post::Vector{Float64}; code = d_code, nucm = CodonMolecularEvolution.demo_nucmat) = dNdS(HB98AA_matrix(1.0, nucm, fs_post, genetic_code = code), HB98AA_matrix(1.0, nucm, zeros(20), genetic_code = code), exp(HB98AA_matrix(1.0, nucm, fs_pre, genetic_code = code) * 100)[1,:], code = code)

"""
    std2maxdNdS(σ)

Approximation for the maximum dN/dS ratio as a function of the standard deviation of the fitnesses, assuming Gaussian fitnesses and a Halpern and Bruno model,
where the fitnesses have just shifted from one Gaussian sample to another. Note: this is not an analytical solution, but a serindipitously good approximation.

```
function monte_carlo_maxdNdS(σ; N=100_000)
    sum_val = 0.0
    for _ in 1:N
        f_i = σ * randn()
        f_j = σ * randn()
        sum_val += HB_fixation_rate(f_i, f_j)
    end
    return sum_val / N
end
vs = 0:0.01:10
plot(vs, monte_carlo_maxdNdS.(vs), label = "Monte Carlo", alpha = 0.8)
plot!(vs, std2maxdNdS.(vs), label = "Approx", linestyle = :dash, alpha = 0.8)
```
"""
std2maxdNdS(σ) = sqrt(σ^2 + π) / sqrt(π)

"""
    maxdNdS2std(ω)

Inverse of std2maxdNdS(σ). Estimates the standard deviation of the fitnesses that will produce, in expectation, a dN/dS ratio of `ω`, assuming Gaussian fitnesses and a Halpern and Bruno model,
where the fitnesses have just shifted from one Gaussian sample to another. Note: this is not an analytical solution, but a serindipitously good approximation.
"""
maxdNdS2std(ω) = sqrt(π * (ω^2 - 1))

"""
    time_varying_HB_freqs(ts, T, fst, init_freqs; nucm = CodonMolecularEvolution.demo_nucmat, alpha = 1.0, delta_t = 0.002, prezero_delta_t = 0.5)

Compute the time-varying codon frequencies and expected dN/dS over time for a sequence of fitnesses, under the Halpern-Bruno model.
`ts` is a vector of times, `T` is the total time, `fst` is a vector of vector of fitnesses, `init_freqs` is the initial codon frequencies,
`nucm` is the nucleotide substitution matrix, `alpha` is the alpha parameter, `delta_t` is the discretization time step for the simulation,
and `prezero_delta_t` is the time step used before `t=0`. `fst[i]` is assumed to be the fitness between `t = ts[i]` and `t = ts[i+1]`.
"""
function time_varying_HB_freqs(ts, T, fst, init_freqs; nucm = CodonMolecularEvolution.demo_nucmat, alpha = 1.0, delta_t = 0.002, prezero_delta_t = 0.5)
    cod_freq_collection = []
    t_collection = []
    dnds_collection = []
    ind = 1
    ts = [ts; Inf]
    t = ts[ind]
    codon_freqs = copy(init_freqs)
    while (t < T) && (ind + 1 <= length(ts))
        if t < 0-prezero_delta_t
            dt = min(prezero_delta_t, ts[ind + 1] - t)
        else
            dt = min(delta_t, ts[ind + 1] - t)
        end
        Q = HB98AA_matrix(alpha, nucm, fst[ind])
        P = exp(Q * dt)
        codon_freqs = (codon_freqs' * P)'
        push!(cod_freq_collection, copy(codon_freqs))
        push!(dnds_collection, dNdS(Q, HB98AA_matrix(alpha, nucm, zeros(20)), codon_freqs))
        if t + dt >= ts[ind + 1]
            ind += 1
        end
        t += dt
        push!(t_collection, t)
    end
    return cod_freq_collection, t_collection, dnds_collection
end


function piecewise_linear_plot!(t_vec, fs_vec, tmax; shift = NaN, colors = nothing, kw...)
    fs = zeros(length(fs_vec[1]), length(fs_vec)*3)
    ts = zeros(length(fs_vec)*3)
    push!(t_vec, tmax)
    for i in 1:length(fs_vec)
        fs[:, 3i-2:3i-1] .= fs_vec[i]
        fs[:, 3i] .= fs_vec[i] .+ shift
        ts[3i-2:3i] .= [t_vec[i], t_vec[i+1], t_vec[i+1]]
    end
    for i in 1:size(fs, 1)
        plot!(ts, fs[i,:], label = :none, color = colors[i]; kw...)
    end
end

function add_muller!(pos_mean_matrix; x_ax = 1:size(pos_mean_matrix)[2], plot_theme = nothing, plot_perm = 1:size(pos_mean_matrix)[1],
    edge_pad = 6, fillalpha = 0.8, colormap = 1:size(pos_mean_matrix)[1])
    pmm = pos_mean_matrix[plot_perm,:]
    cum_freq = zeros(size(pmm)[2])
    for i in 1:size(pmm)[1]
        plot!(x_ax,1 .- cum_freq, fillrange = 1 .- (cum_freq .+ pmm[i,:]), fillalpha = fillalpha, linewidth = 0, color = plot_theme[colormap[plot_perm[i]]], label = :none)
        cum_freq .+= pmm[i,:]
    end
end

"""
    HBviz(ts, fst, T, alp, nucm)

Visualize over time the fitness trajectory, the codon frequencies, and the expected dN/dS. `ts` is a vector of times, `fst` is a vector of fitnesses, `T` is the total time, `alp` is the alpha parameter, and `nucm` is the nucleotide substitution matrix.

```julia
σ = 2.0
alpha = 1.0
nucm = CodonMolecularEvolution.demo_nucmat
fst = [randn(20) .* σ, randn(20) .* σ]
ts = [-100.0, 1.0]
T = 2.0
HBviz(ts, fst, T, alpha, nucm)
```
"""
function HBviz(ts::Vector{Float64}, fst::Vector{Vector{Float64}}, T::Float64, alp, nucm;
                scram = [8, 4, 11, 17, 9, 14, 18, 2, 10, 13, 16, 20, 12, 6, 19, 7, 1, 15, 5, 3],
                viz_size = (1200,500), σ = nothing)
    cfc, tc, dndses = time_varying_HB_freqs(ts, T, fst, ones(61)/61, alpha = alp)
    #Removing burnin:
    prezero = findlast(tc .<= 0)
    prezero = prezero == nothing ? 1 : prezero
    tc = tc[prezero:end]
    cfc = cfc[prezero:end]
    dndses = dndses[prezero:end]
    fl = findlast(ts .<= 0)
    fl = fl == nothing ? 1 : fl
    ts = ts[fl:end]
    fst = fst[fl:end]
    ###
    perm = sortperm(MolecularEvolution.universal_code.amino_acid_lookup)
    AA_nums = MolecularEvolution.universal_code.codon2AA_pos[perm]
    colorv = zeros(Int, 61)
    c = 1
    for i in 1:20
        s = findall(AA_nums .== i)
        colorv[c:c+length(s)-1] .= 10*(scram[i]-1) .+ collect(1:length(s)) .* 2 .- 1
        c += length(s)
    end
    plot_theme_exploded = cgrad(:darkrainbow, 200, categorical = true)
    AA_plot_theme = plot_theme_exploded[1 .+ (10 .* (scram .- 1))]
    plot_theme = plot_theme_exploded[colorv]
    pl1 = plot(tc, zeros(length(tc)), size = viz_size, xlim = (0, T), linestyle = :dash, color = "black", alpha = 0.0,
                ylabel = L"f_t", xtickfontcolor = RGBA(0,0,0,0), legend=false)
    extr = maximum([maximum(abs.(f)) for f in fst]) * 1.1
    piecewise_linear_plot!(ts, fst, T, colors = AA_plot_theme, ylims = (-extr, extr))
    pl2 = plot(tc, zeros(length(tc)), size = viz_size, alpha = 0.0, xlim = (0, T), ylabel = L"C_t", xtickfontcolor = RGBA(0,0,0,0), legend=false, top_margin = -12Plots.mm)
    add_muller!(stack(cfc)[perm,:], plot_theme = plot_theme, x_ax = tc)
    pl3 = plot(tc, dndses, size = viz_size, xlim = (0, T), xlabel = "Time", ylabel = L"E_C_t[dN/dS|f_t]", ylim = (0, maximum(dndses[tc .>= 0]) + 0.1),
        label = :none, top_margin = -12Plots.mm, color = "black")
    plot!([0,T], [1,1], color = "black", linestyle = :dash, alpha = 0.5, label = :none)
    if !isnothing(σ)
        maxdnds = std2maxdNdS(σ)
        plot!([0,T], [maxdnds,maxdnds], color = "red", linestyle = :dash, alpha = 0.5, label = :none)
    end
    return plot(pl1, pl2, pl3, layout=(3, 1), link=:x, margins = 8Plots.mm, plot_layout = :tight, widen=false, tickdirection=:out)
end

"""
    shiftingHBviz(T, event_rate, σ, mixing_rate, alpha, nucm; T0 = -20)

Visualize the fitness trajectory, codon frequencies, and expected dN/dS over time for a shifting HB process.
`T` is the total time, `event_rate` is the rate of fitness shifts, `σ` is the standard deviation of the fitnesses,
`mixing_rate` is the rate of mixing between fitnesses, `alpha` is the alpha parameter, and `nucm` is the nucleotide substitution matrix.
`T0` controls the burnin time, to ensure the process is at equilibrium at `t=0`.

```julia
T = 2.0
mix = 1.0
σ = 5.0
event_rate = 100.0
alpha = 1.0
nucm = CodonMolecularEvolution.demo_nucmat
shiftingHBviz(T, event_rate, σ, mix, alpha, nucm)
```
"""
function shiftingHBviz(T, event_rate, σ, mixing_rate, alpha, nucm; T0 = -20)
    fs = randn(20) .* σ
    m = PiecewiseOUModel(event_rate, σ, mixing_rate)
    coll = []
    codon = sample(1:61, Weights(HB98AA_eqfreqs(fs)))
    newfs, _, _, jumps = jumpy_HB_codon_evolve(fs, codon, m, nucm, alpha, T-T0, push_into = coll)
    prepend!(coll, [(0.0, codon, fs)])
    ts = [c[1] for c in coll] .+ T0
    fst = [c[3] for c in coll]
    return HBviz(ts, fst, T, alpha, nucm, σ = σ)
end



###############################################################################
#A model for evolving log effective pop size & unscaled site fitnesses
###############################################################################

# jumpy_NeHB_codon_evolve is currently fine for a single site. However, to get this version working over a
#branch for multiple sites, we need to:
# 1) First sample the Ne jumps over the branch (since these apply for all sites)
# 2) One site at a time, sample the fitness jumps and codon substitutions between each Ne jump event
# ShiftingNeHBSimModel and ShiftingNeHBSimPartition are commented out until this works.
"""
    jumpy_NeHB_codon_evolve(fitnesses, logNe_trajectory, codon, fitness_model, nuc_matrix, alpha;
    genetic_code = MolecularEvolution.universal_code, push_into = nothing)

Evolves codons and unscaled site-fitness, along with a given trajectory of log-pop-size.
"""
function jumpy_NeHB_codon_evolve(fitnesses, logNe_trajectory, codon, fitness_model, nuc_matrix, alpha;
    genetic_code = MolecularEvolution.universal_code, push_into = nothing)
    codon_jumps = 0
    fitness_jumps = 0
    current_fits = copy(fitnesses)
    current_codon = codon
    t = 0.0
    for (next_event, current_logNe) in logNe_trajectory
        current_fits, current_codon, local_codon_jumps, local_fitness_jumps = jumpy_HB_codon_evolve(current_fits, current_codon, fitness_model, nuc_matrix, alpha, next_event; genetic_code = genetic_code, push_into = push_into, logNe = current_logNe, time_origin = t)
        t = t+next_event
        codon_jumps += local_codon_jumps
        fitness_jumps += local_fitness_jumps
        if !isnothing(push_into)
            push!(push_into,(t,current_codon,copy(current_fits) .* exp(current_logNe) .* 2,current_logNe))
        end
    end
    return logNe_trajectory[end][2], current_fits, current_codon, codon_jumps, fitness_jumps, length(logNe_trajectory)-1
end


function jumpy_Ne(logNe, logNe_model, time)
    current_logNe = logNe
    t = 0.0
    next_event = 0.0
    logNe_trajectory = Vector{Tuple{Float64, Float64}}()
    while t+next_event < time
        next_event = randexp()/logNe_model.event_rate
        t = t+next_event
        if t < time
            push!(logNe_trajectory, (next_event, current_logNe))
            current_logNe = jump(current_logNe, logNe_model)
        end
    end
    push!(logNe_trajectory, (time - (t - next_event), current_logNe))
    return logNe_trajectory
end

"""
    shiftingNeHBviz(T, f_event_rate, f_σ, f_mixing_rate, logNe_event_rate, logNe_σ, logNe_mean, logNe_mixing_rate, alpha, nucm; T0 = -20)

Visualize the Ne trajectory, fitness trajectory, codon frequencies, and expected dN/dS over time for a shifting Ne HB process.
"""
function shiftingNeHBviz(T, f_event_rate, f_σ, f_mixing_rate, logNe_event_rate, logNe_σ, logNe_mean, logNe_mixing_rate, alpha, nucm; T0 = -20)
    fs = randn(20) .* f_σ
    logNe = randn() .* logNe_σ .+ logNe_mean
    f_ou = PiecewiseOUModel(f_event_rate, f_σ, f_mixing_rate)
    log_ne_ou = PiecewiseOUModel(logNe_event_rate, logNe_σ, logNe_mixing_rate, logNe_mean, 0.0, 0.0)
    coll = []
    codon = sample(1:61, Weights(HB98AA_eqfreqs(fs .* exp(logNe) .* 2)))
    logNe_trajectory = jumpy_Ne(logNe, log_ne_ou, T-T0)
    jumpy_NeHB_codon_evolve(fs, logNe_trajectory, codon, f_ou, nucm, alpha, push_into = coll)
    prepend!(coll, [(0.0, codon, fs, logNe)])
    ts = [c[1] for c in coll] .+ T0
    fst = [c[3] .* exp(c[4]) .* 2 for c in coll] #2Ne*s
    p = HBviz(ts, fst, T, alpha, nucm)

    prezero = findlast(ts .<= 0)
    prezero = prezero == nothing ? 1 : prezero
    pl4 = plot(ts[prezero:end], exp.([c[4] for c in coll[prezero:end]]), ylabel = L"Ne", legend = :none, color = "black", xtickfontcolor = RGBA(0,0,0,0), bottom_margin = -12Plots.mm, xlim = (0, T))
    return plot(pl4, p, layout = grid(2, 1, heights = 1/4 .*[1,3]), link=:x, margins = 8Plots.mm, plot_layout = :tight, widen=false, tickdirection=:out)
end


mutable struct ShiftingNeHBSimModel <: MolecularEvolution.SimulationModel
    sites::Int64
    alphas::Vector{Float64}                # per-site synonymous rates
    unscaled_ou_params::Vector{PiecewiseOUModel}  # OU process for each site's unscaled fitness
    logNe_model::PiecewiseOUModel       # OU process for log(pop size)
    nuc_matrix::Matrix{Float64}
    code::MolecularEvolution.GeneticCode
end

"""
    ShiftingNeHBSimModel(nuc_matrix, unscaled_ou_params, logNe_model; alpha, rescale, code)

Create a model with one OU process per site (unscaled_ou_params) and one OU process for log(N_e). 
`alpha` can be a scalar or a vector of the same length as unscaled_ou_params.
If `rescale` is true, the nuc_matrix is scaled (HB98 neutral scaling).
"""
function ShiftingNeHBSimModel(
    nuc_matrix::Matrix{Float64}, 
    unscaled_ou_params::Vector{PiecewiseOUModel{A,B,C}}, 
    logNe_model::PiecewiseOUModel;
    alpha = ones(length(unscaled_ou_params)),
    rescale::Bool = true,
    code::MolecularEvolution.GeneticCode = MolecularEvolution.universal_code) where {A,B,C}
    sites = length(unscaled_ou_params)
    @assert length(alpha) == sites || length(alpha) == 1 "alpha must match number of sites or be a single value"

    # Copy or scale the nuc matrix if requested
    mat = rescale ? HB98AA_neutral_scale_nuc_matrix(copy(nuc_matrix), genetic_code = code) : copy(nuc_matrix)

    return ShiftingNeHBSimModel(
        sites,
        length(alpha) == 1 ? fill(alpha[1], sites) : alpha,
        unscaled_ou_params,
        logNe_model,
        mat,
        code
    )
end

mutable struct ShiftingNeHBSimPartition <: CodonSimulationPartition
    sites::Int64
    logNe::Float64
    fitnesses::Matrix{Float64}  # size: (number_of_AAs, sites)
    codons::Vector{Int64}
    code::MolecularEvolution.GeneticCode
end

"""
    ShiftingNeHBSimPartition(model; burnin_time)

Constructs a partition for the `ShiftingNeHBSimModel`.
Initializes by:
  - Drawing random unscaled fitnesses from each site's OU equilibrium.
  - Drawing random logNe from its OU equilibrium.
  - "Burning in" each site along the branch of length `burnin_time`.
"""
function ShiftingNeHBSimPartition(
    model::ShiftingNeHBSimModel; 
    burnin_time::Float64 = 100.0
)
    # 1) sample an initial logNe from equilibrium
    # The equilibrium of an OU in this piecewise-constant approximation is roughly Normal(mean=mu, sd=eq_std).
    lp = randn() * model.logNe_model.eq_std + model.logNe_model.mu

    # 2) sample initial unscaled fitness for each site
    nAAs = length(model.code.amino_acids)
    fits = zeros(nAAs, model.sites)
    for s in 1:model.sites
        oup = model.unscaled_ou_params[s]
        fits[:,s] .= randn(nAAs) .* oup.eq_std .+ oup.mu
    end

    # 3) pick initial codons (arbitrary, say the first sense codon)
    codons = [rand(1:length(model.code.sense_codons)) for i in 1:model.sites]

    # 4) "burn in" by evolving over burnin_time
    #jumpy_NeHB_codon_evolve(fitnesses, logNe, codon, fitness_model, logNe_model, nuc_matrix, alpha, time
    logNe_trajectory = jumpy_Ne(lp, model.logNe_model, burnin_time)
    for (i,m) in enumerate(model.unscaled_ou_params)
        _, f, c, _, _ = jumpy_NeHB_codon_evolve(fits[:,i], logNe_trajectory, codons[i], m, model.nuc_matrix, model.alphas[i], genetic_code = model.code)
        fits[:,i] .= f
        codons[i] = c
    end
    logNe = logNe_trajectory[end][2]

    return ShiftingNeHBSimPartition(
        model.sites,
        logNe,
        fits,
        codons,
        model.code
    )
end

# Minimal constructor for placeholders (if one needs a zero-time partition, for instance)
function ShiftingNeHBSimPartition(
    sites::Int; 
    code::MolecularEvolution.GeneticCode = MolecularEvolution.universal_code
)
    nAAs = length(code.amino_acids)
    return ShiftingNeHBSimPartition(
        sites,
        0.0, 
        zeros(nAAs, sites),
        ones(Int, sites),
        code
    )
end

"""
    forward!(dest, source, model, node)

Evolves `source` partition along `node.branchlength` under `model`, storing the result in `dest`.
"""
function MolecularEvolution.forward!(
    dest::ShiftingNeHBSimPartition,
    source::ShiftingNeHBSimPartition,
    model::ShiftingNeHBSimModel,
    node::FelNode
)
    bl = node.branchlength
    logNe_trajectory = jumpy_Ne(source.logNe, model.logNe_model, bl)
    logNe = logNe_trajectory[end][2]
    for site in 1:model.sites
        fitnesses = source.fitnesses[:,site]
        codon = source.codons[site]
        unscaled_ou_model = model.unscaled_ou_params[site]
        alpha = model.alphas[site]
        _, fitnesses, codon, _, _, _ = jumpy_NeHB_codon_evolve(
            fitnesses,
            logNe_trajectory,
            codon,
            unscaled_ou_model,
            model.nuc_matrix,
            alpha;
            genetic_code = model.code
        )
        dest.fitnesses[:,site] .= fitnesses
        dest.codons[site] = codon
    end
    dest.logNe = logNe
end

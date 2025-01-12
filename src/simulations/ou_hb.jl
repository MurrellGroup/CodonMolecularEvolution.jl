#Halpern and Bruno model where amino acid fitnesses evolve over time using a piecewise constant approximation to an OU process.
#Authors: Hassan Sadiq and Ben Murrell

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
    PiecewiseOUModel(event_rate::Float64, eq_std::Float64,independence::Float64; delta_t = 1.0)
    PiecewiseOUModel(offsets::Vector{Float64})

A piecewise constant approximation to an OU process, intended to simulate fitnesses evolving over phylogenies.
The equilibrium standard deviation is directly parameterized (`eq_std`), as is the rate at which the process mixes to equilibrium (`mixing`).
`event_rate` controls how often the fitness changes occur, where the mixing rate is scaled to compensate for the increased rate of change to achieve
approximately the same amount of change per unit time even as the `event_rate` changes. A very high `event_rate` will behave more like continuous diffusion,
but will be more computationally expensive to sample from. `mu` can also be set to control the mean fitnesses.
The model also permits `offsets`, which are added to the fitnesses as they are passed into the model. For a single process, these are confounded with the mean `mu`
but if the offsets change (eg. from one branch to another) the effective fitnesses will immidiately change, whereas if `mu` changes the fitnesses will drift towards `mu`.
"""
mutable struct PiecewiseOUModel{M<:Union{Float64,AbstractVector{Float64}},O<:Union{Float64,AbstractVector{Float64}}} 
    mu::M #OU mean. Scalar (in which case it is the same for all AAs) or vector (one per AA)
    delta_t::Float64 #Scales the size of the fitness jumps
    var::Float64 #OU variance
    reversion_rate::Float64 #Rate of reversion to the mean
    event_rate::Float64 #Rate of fitness jumps
    offsets::O #These do not change, thus allowing preferred codons. Incorporated right at the end.
    eq_std::Float64 #Standard deviation of the equilibrium distribution
    function PiecewiseOUModel(event_rate::Float64, eq_std::Float64, mixing::Float64; delta_t = 1.0)
        new{Float64,Float64}(0.0, 1.0, 2*(mixing/event_rate)*(eq_std^2), mixing/event_rate, event_rate, 0.0, eq_std)
    end
    function PiecewiseOUModel(offsets::Vector{Float64}) #Constructor for static model, using just offsets for fitnesses
        new{Float64,Vector{Float64}}(0.0, 1.0, 1e-15, 1.0, 0.0, offsets, 1e-15)
    end
end

"""
    ou_jump(x_source, delta_t::Float64, mu::Union{Float64,Vector{Float64}}, sigma2::Float64, theta::Float64)

Evolves values over time using a piecewise constant approximation to an OU process, where this function computes the new distribution for a single discrete jump.
`x_source` is the vector of fitnesses, `delta_t` is the time step, `mu` is the mean fitness, `sigma2` is the variance, and `theta` is the reversion rate.
"""
function ou_jump(x_source, delta_t::Float64, mu::Union{Float64,Vector{Float64}}, sigma2::Float64, theta::Float64)
    sqrt_ou_var = sqrt((sigma2/(2*theta)) * (1 - exp(-2 * theta * delta_t)))
    mu_multiplier = exp(-theta * delta_t)
    return randn(length(x_source)) .* sqrt_ou_var .+ mu .+ ((x_source .- mu) .* mu_multiplier)
end

function ou_jump(x_source, m::PiecewiseOUModel)
    return ou_jump(x_source, m.delta_t, m.mu, m.var, m.reversion_rate)
end

#Todo: Introduce a model with jumps to random equilibrium fits using dispatch here, just by overloading ou_jump on a different model type.

"""
    HB98AA_row(current_codon, alpha, nuc_matrix, AA_fitness; genetic_code=MolecularEvolution.universal_code)

Returns the rate row for a codon model using the HB98 model where each AA has a different fitness. `current_codon` is the current codon, `alpha` is the synonymous rate,
`nuc_matrix` is the symmetric nucleotide substitution matrix, and `AA_fitness` is the fitness of each amino acid.
"""
function HB98AA_row(current_codon, alpha, nuc_matrix, AA_fitness; genetic_code=MolecularEvolution.universal_code)
    row = zeros(genetic_code.num_sense)
    codon_aa_i = AA_fitness[genetic_code.codon2AA_pos[current_codon]]
    for j in 1:genetic_code.num_sense
        if genetic_code.codon_to_nuc_map[current_codon,j] != (-1,-1,-1)
            c2n_map = genetic_code.codon_to_nuc_map[current_codon,j]
            f_ab = HB_fixation_rate(codon_aa_i,AA_fitness[genetic_code.codon2AA_pos[j]])
            row[j] = alpha * nuc_matrix[c2n_map[2],c2n_map[3]] * f_ab
        end
    end
    return row
end

"""
    jumpy_HB_codon_evolve(fitnesses, codon, ou_model, nuc_matrix, alpha, time;
        genetic_code = MolecularEvolution.universal_code, push_into = nothing)

Evolves fitnesses and codons over time using the HB98 model. `fitnesses` is the vector of fitnesses, `codon` is the current codon, `ou_model` is the OU model,
`nuc_matrix` is the symmetric nucleotide substitution matrix, `alpha` is the synonymous rate, and `time` is the total time to evolve over.
"""
function jumpy_HB_codon_evolve(fitnesses, codon, ou_model, nuc_matrix, alpha, time;
        genetic_code = MolecularEvolution.universal_code, push_into = nothing)
    codon_jumps = 0
    fitness_jumps = 0
    current_fits = copy(fitnesses)
    current_codon = codon
    t = 0.0
    next_event = 0.0
    while t+next_event < time
        HBrow = HB98AA_row(current_codon, alpha, nuc_matrix, current_fits .+ ou_model.offsets, genetic_code=genetic_code)
        sum_HBrow = sum(HBrow)
        rOU,rHB = (ou_model.event_rate,sum_HBrow)
        total_rate = rOU+rHB
        next_event = randexp()/total_rate
        t = t+next_event
        if t < time
            event_index = sample(1:2,Weights([rOU,rHB])) 
            if event_index == 1 # Fitness jump event
                fitness_jumps += 1
                current_fits = ou_jump(current_fits, ou_model)
            else # Codon substitution event
                codon_jumps += 1
                current_codon = sample(1:length(HBrow),Weights(HBrow))
            end
        end
        if !isnothing(push_into)
            push!(push_into,(t,current_codon,copy(current_fits)))
        end
    end
    return current_fits, current_codon, codon_jumps, fitness_jumps
end


#=
function piecewise_linear_plot!(t_vec, fs_vec, tmax; shift = NaN, kw...)
    fs = zeros(length(fs_vec[1]), length(fs_vec)*3)
    ts = zeros(length(fs_vec)*3)
    push!(t_vec, tmax)
    for i in 1:length(fs_vec)
        fs[:, 3i-2:3i-1] .= fs_vec[i]
        fs[:, 3i] .= fs_vec[i] .+ shift
        ts[3i-2:3i] .= [t_vec[i], t_vec[i+1], t_vec[i+1]]
    end
    for i in 1:size(fs, 1)
        plot!(ts, fs[i,:], label = :none; kw...)
    end
end

#For investigating mixing rates:
t = 1.0 #Try 1.0 and 20.0
fs = zeros(20) .+ 5
pl = plot(xlim = (0,t))
rates = [0.2, 1.0, 5.0, 25.0]
colors = ["red", "blue", "green", "purple"]
coll = nothing
for i in 1:4
    for _ in 1:4
        m = PiecewiseOUModel(1000.0/t, 1.0, rates[i])
        coll = []
        newfs, _, _, jumps = jumpy_HB_codon_evolve(fs, 1, m, CodonMolecularEvolution.demo_nucmat, 0.0, t, push_into = coll)
        prepend!(coll, [(0.0, 1, fs)])
        piecewise_linear_plot!([c[1] for c in coll], [c[3] for c in coll], t, color = colors[i], alpha = 0.1)
    end
    piecewise_linear_plot!([c[1] for c in coll], [c[3][1:1] for c in coll], t, color = colors[i], alpha = 0.7, shift = 0.0, label = "$(rates[i])")
end
pl

#Confirming that the jump distribution compensates for increasing the event rate:
t = 1.0
fs = zeros(20) .+ 5
pl = plot(xlim = (0,t))
m = PiecewiseOUModel(1000.0, 1.0, 5.0)
for i in 1:4
    coll = []
    newfs, _, _, jumps = jumpy_HB_codon_evolve(fs, 1, m, CodonMolecularEvolution.demo_nucmat, 0.0, t, push_into = coll)
    prepend!(coll, [(0.0, 1, fs)])
    piecewise_linear_plot!([c[1] for c in coll], [c[3] for c in coll], t, color = "red", alpha = 0.1)
end
m = PiecewiseOUModel(100.0, 1.0, 5.0)
for i in 1:4
    coll = []
    newfs, _, _, jumps = jumpy_HB_codon_evolve(fs, 1, m, CodonMolecularEvolution.demo_nucmat, 0.0, t, push_into = coll)
    prepend!(coll, [(0.0, 1, fs)])
    piecewise_linear_plot!([c[1] for c in coll], [c[3] for c in coll], t, color = "blue", alpha = 0.1)
end
pl
=#

"""
    ShiftingHBSimPartition(nuc_matrix::Matrix{Float64}, models::Vector{PiecewiseOUModel}; burnin_time = 100.0, code = MolecularEvolution.universal_code)

Constructs a partition that tracks evolving fitnesses and codons. Only useable for sampling (not likelihood calculations).
"""
mutable struct ShiftingHBSimPartition <: Partition
    sites::Int64
    fitnesses::Matrix{Float64}
    codons::Vector{Int64}
    code::MolecularEvolution.GeneticCode
    function ShiftingHBSimPartition(sites; code = MolecularEvolution.universal_code)
        new(sites,zeros(length(code.amino_acids),sites),ones(Int64,sites),code)
    end
    function ShiftingHBSimPartition(nuc_matrix::Matrix{Float64}, models::Vector{PiecewiseOUModel}; burnin_time = 100.0, code = MolecularEvolution.universal_code)
        fits = zeros(length(code.amino_acids),length(models))
        for (i,m) in enumerate(models)
            fits[:,i] .= randn(length(code.amino_acids)) .* m.eq_std .+ m.mu
        end
        codons = ones(Int64,length(models))
        for (i,m) in enumerate(models)
            f, c, _, _ = jumpy_HB_codon_evolve(fits[:,i], codons[i], m, nuc_matrix, 1.0, burnin_time, genetic_code = code)
            fits[:,i] .= f
            codons[i] = c
        end
        new(length(models), fits, codons)
    end
end

"""
    ShiftingHBSimModel(sites, alphas, ou_params, nuc_matrix)

A model for simulating fitnesses evolving over phylogenies using the HB98 model. `sites` is the number of sites, `alphas` is a vector of synonymous rates (one per site),
`ou_params` is a vector of `PiecewiseOUModel`s (one per site), and `nuc_matrix` is the symmetric nucleotide substitution matrix (shared across sites).
"""
mutable struct ShiftingHBSimModel <: MolecularEvolution.SimulationModel
    sites::Int64
    alphas::Vector{Float64}
    ou_params::Vector{PiecewiseOUModel}
    nuc_matrix::Matrix{Float64}
end

#Does nothing because the `forward!` function implicitly samples.
function MolecularEvolution.sample_partition!(p::ShiftingHBSimPartition)
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
            node.branchlength)
        dest.fitnesses[:,site] .= fitnesses
        dest.codons[site] = codon
    end
end

function partition2string(part::ShiftingHBSimPartition; code = MolecularEvolution.universal_code)
    return join([code.sense_codons[part.codons[i]] for i in 1:part.sites])
end

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

#=
fs = randn(20)
c = rand(1:61)
r1 = HB98AA_row(c, 1.23, CodonMolecularEvolution.demo_nucmat, fs)
q1 = HB98AA_matrix(1.23, CodonMolecularEvolution.demo_nucmat, fs)
q1[c,c] = 0
@assert isapprox(sum(q1[c,:] .- r1), 0.0)
=#

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
function dNdS(q1, q0, p; code = d_code)
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
    dNdS(q1, q0, exp(q1 * 100)[1,:], code = code)
end

HBdNdS(fs_pre::Vector{Float64}, fs_post::Vector{Float64}; code = d_code, nucm = CodonMolecularEvolution.demo_nucmat) = dNdS(HB98AA_matrix(1.0, nucm, fs_post, genetic_code = code), HB98AA_matrix(1.0, nucm, zeros(20), genetic_code = code), exp(HB98AA_matrix(1.0, nucm, fs_pre, genetic_code = code) * 100)[1,:], code = code)

"""
    approx_std2maxdNdS(σ)

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
plot!(vs, approx_std2maxdNdS.(vs), label = "Approx", linestyle = :dash, alpha = 0.8)
```
"""
approx_std2maxdNdS(σ) = sqrt(σ^2 + π) / sqrt(π)

"""
    approx_maxdNdS2std(ω)

Inverse of approx_std2maxdNdS(σ). Estimates the standard deviation of the fitnesses that will produce, in expectation, a dN/dS ratio of `ω`, assuming Gaussian fitnesses and a Halpern and Bruno model,
where the fitnesses have just shifted from one Gaussian sample to another. Note: this is not an analytical solution, but a serindipitously good approximation.
"""
approx_maxdNdS2std(ω) = sqrt(π * (ω^2 - 1))



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
    perm = sortperm(MolecularEvolution.universal_code.amino_acid_lookup)
    AA_nums = MolecularEvolution.universal_code.codon2AA_pos[perm]
    colorv = zeros(Int, 61)
    c = 1
    for i in 1:20
        s = findall(AA_nums .== i)
        colorv[c:c+length(s)-1] .= 10*(scram[i]-1) .+ collect(1:length(s)) .* 2 .- 1
        c += length(s)
    end
    plot_theme_exploded = cgrad(:rainbow, 200, categorical = true)
    AA_plot_theme = plot_theme_exploded[1 .+ (10 .* (scram .- 1))]
    plot_theme = plot_theme_exploded[colorv]
    pl1 = plot(tc, zeros(length(tc)), size = viz_size, xlim = (0, T), linestyle = :dash, color = "black", alpha = 0.0,
                ylabel = L"f_t", xtickfontcolor = RGBA(0,0,0,0), legend=false)
    fl = findlast(ts .<= 0)
    extr = maximum([maximum(abs.(f)) for f in fst[fl:end]]) * 1.1
    piecewise_linear_plot!(ts, fst, T, colors = AA_plot_theme, ylims = (-extr, extr))
    pl2 = plot(tc, zeros(length(tc)), size = viz_size, alpha = 0.0, xlim = (0, T), ylabel = L"C_t", xtickfontcolor = RGBA(0,0,0,0), legend=false, top_margin = -12Plots.mm)
    add_muller!(stack(cfc)[perm,:], plot_theme = plot_theme, x_ax = tc)
    pl3 = plot(tc, dndses, size = viz_size, xlim = (0, T), xlabel = "Time", ylabel = L"E_C_t[dN/dS|f_t]", ylim = (0, maximum(dndses[tc .>= 0]) + 0.1),
        label = :none, top_margin = -12Plots.mm, color = "black")
    plot!([0,T], [1,1], color = "black", linestyle = :dash, alpha = 0.5, label = :none)
    if !isnothing(σ)
        maxdnds = approx_std2maxdNdS(σ)
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
    newfs, _, _, jumps = jumpy_HB_codon_evolve(fs, 1, m, nucm, alpha, T-T0, push_into = coll)
    prepend!(coll, [(0.0, 1, fs)])
    ts = [c[1] for c in coll] .+ T0
    fst = [c[3] for c in coll]
    return HBviz(ts, fst, T, alpha, nucm, σ = σ)
end

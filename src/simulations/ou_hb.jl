#Halpern and Bruno model where amino acid fitnesses evolve over time using a piecewise constant approximation to an OU process.
#Authors: Hassan Sadiq and Ben Murrell

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
    dNdS(q1, q0; code = MolecularEvolution.universal_code)

Returns an analytic estimate of the dN/dS ratio for a codon model. `q1` is the rate matrix where selection is active (eg. a Halpern and Bruno model with a set of fitnesses),
and `q0` is the corresponding rate matrix where selection is inactive (eg. a Halpern and Bruno model with all fitnesses equal).

```
fs = randn(20)
nucm = CodonMolecularEvolution.demo_nucmat
q1 = HB98AA_matrix(1.0, nucm, fs)
q0 = HB98AA_matrix(1.0, nucm, zeros(20))
dNdS(q1, q0)
```
"""
function dNdS(q1, q0; code = MolecularEvolution.universal_code)
    p = exp(q1 * 100)[1,:]
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

function MEME_global_fit(seqnames::Vector{String}, seqs, treestring::String;
    verbosity=1, code=MolecularEvolution.universal_code, optimize_branch_lengths=false)
    tree = FUBAR_init(treestring, verbosity=verbosity)
    difFUBAR_global_fit_2steps(seqnames, seqs, tree, x -> x, code, verbosity=verbosity, optimize_branch_lengths=optimize_branch_lengths)
end

function build_model_vec(p, GTRmat, F3x4, genetic_code)
    return [GeneralCTMC(MolecularEvolution.MG94_F3x4(p.alpha, p.alpha*o, GTRmat, F3x4, genetic_code=genetic_code)) for o in p.omegas]
end

function build_mixture_model(p, GTRmat, F3x4, genetic_code)
    models = build_model_vec(p, GTRmat, F3x4, genetic_code)
    weights = [p.qminus, 1-p.qminus]
    return BWMModel{eltype(models)}(models, weights)
end

function MEME_site_model_fit(site::Int, omega_restrictions::NamedTuple, initial_params::NamedTuple, tree::FelNode, first_site::Bool, seqnames, seqs, GTRmat, F3x4_freqs, eq_freqs, genetic_code)
    #Fit model for site
    #This should work for both the alternative and null, governed by omega_restriction
    seq_ind = 1 + 3 * (site-1)
    parent_part = tree.parent_message[1]
    if first_site
        eq_partition = CodonPartition(1)
        eq_partition.state .= eq_freqs
        initial_partition = LazyPartition{CodonPartition}()
    else
        initial_partition = parent_part
    end
    populate_tree!(tree, initial_partition, seqnames, [x[seq_ind:seq_ind+2] for x in seqs], init_all_messages = first_site)
    if first_site
        lazyprep!(tree, eq_partition)
    else
        push!(parent_part.memoryblocks, parent_part.partition) #Don't want to renew parent_message every time..., could just pop it hehe, but not very clean...
    end
    tree.parent_message[1].static = true
    
    flat_initial_params, unflatten = value_flatten(initial_params) #See ParameterHandling.jl docs
    num_params = length(flat_initial_params)

    function objective(params::NamedTuple; tree=tree, eq_freqs=eq_freqs)
        return -log_likelihood!(tree, build_mixture_model(params, GTRmat, F3x4_freqs, genetic_code))
    end

    #BOBYQA? we'll try
    opt = Opt(:LN_BOBYQA, num_params)
    min_objective!(opt, (x, y) -> (objective ∘ unflatten)(x))
    lower_bounds!(opt, vcat(omega_restrictions.lower_bounds, 0.0, -5.0))
    upper_bounds!(opt, vcat(omega_restrictions.upper_bounds, 1.0, 5.0))
    xtol_rel!(opt, 1e-12)
    _, mini, _ = NLopt.optimize(opt, flat_initial_params)

    final_params = unflatten(mini)
    final_LL = -objective(final_params)
    return final_params, final_LL
end

function MEME_hypotheses_fit(tree::FelNode, seqnames, seqs, GTRmat, F3x4_freqs, eq_freqs, genetic_code; verbosity = 1)
    #Fit alternative model and null model for every site
    num_sites = sites(tree.parent_message[1])
    alternative_params = Vector{NamedTuple}(undef, num_sites)
    null_params = Vector{NamedTuple}(undef, num_sites)
    alternative_LLs = Vector{Float64}(undef, num_sites)
    null_LLs = Vector{Float64}(undef, num_sites)
    #Bounds are in log-domain, except qminus
    alternative_omega_restrictions = (lower_bounds = [-10.0, -10.0], upper_bounds = [0.0, 10.0]) #exp(10) ≈ unrestricted
    null_omega_restrictions = (lower_bounds = [-10.0, -10.0], upper_bounds = [0.0, 0.0])
    global_initial_params = (
        omegas=positive([0.5, 2.0]),
        qminus=0.5,
        alpha=positive(1.0)
    )
    iter = 1:num_sites
    for site in iter #Can easily be parallelized over sites
        #Can make this cleaner with a struct or smth
        for (isnull, params_vec, LLs_vec, omega_restrictions) in zip([false, true], [alternative_params, null_params], [alternative_LLs, null_LLs], [alternative_omega_restrictions, null_omega_restrictions])
            if isnull
                final_alt_params = alternative_params[site]
                initial_params = (
                    omegas=positive([final_alt_params.omegas[1], min(1.0, final_alt_params.omegas[2])]), #Use old optim if alt beta+ < alpha
                    qminus=final_alt_params.qminus,
                    alpha=positive(final_alt_params.alpha)
                )
            else
                initial_params = global_initial_params
            end
            final_params, final_LL = MEME_site_model_fit(
                                        site,
                                        omega_restrictions,
                                        initial_params,
                                        tree,
                                        site == first(iter),
                                        seqnames,
                                        seqs,
                                        GTRmat,
                                        F3x4_freqs,
                                        eq_freqs,
                                        genetic_code
                                     )
            params_vec[site] = final_params
            LLs_vec[site] = final_LL
        end
        verbosity > 0 && if mod(site,5)==1
            print(round(100*site/num_sites),"% ")
            flush(stdout)
        end
    end
    return alternative_params, alternative_LLs, null_params, null_LLs
end

p_value(λ) = 1 - (0.33 + 0.3*cdf(Chisq(1), λ) + 0.37*cdf(Chisq(2), λ)) #Asymptotic (null) distribution of LRT

function MEME_test(alternative_LLs, null_LLs; significance=0.05)
    num_sites = length(alternative_LLs)
    LRTs = Vector{Float64}(undef, num_sites)
    p_values = Vector{Float64}(undef, num_sites)
    for (site, (alternative_LL, null_LL)) in enumerate(zip(alternative_LLs, null_LLs))
        λ = 2 * (alternative_LL - null_LL)
        p = p_value(λ)
        LRTs[site] = λ
        p_values[site] = p
        if p < significance
            #Report that some branches have episodic diversifying selection at site
            println("Site $(site): p-value of the LRT = $(round(p,digits=4))");
        end
    end
    return LRTs, p_values
end

function MEME_tabulate(LRTs, p_values, alternative_params, alternative_LLs, outpath, exports)
    df = DataFrame()
    df."site" = 1:length(p_values)
    df."α" = [alt.alpha for alt in alternative_params]
    df."β⁻" = [alt.alpha * alt.omegas[1] for alt in alternative_params]
    df."q⁻" = [alt.qminus for alt in alternative_params]
    df."β⁺" = [alt.alpha * alt.omegas[2] for alt in alternative_params]
    df."q⁺" = [1 - alt.qminus for alt in alternative_params]
    df."LRT" = LRTs
    df."p-value" = p_values
    df."MEME LogL" = alternative_LLs

    # Round all numerical values to 2 decimal places
    for col in names(df)
        if eltype(df[!, col]) <: AbstractFloat
            df[!, col] = round.(df[!, col], digits=2)
        end
    end
    exports && CSV.write(outpath*"_SelectionOutput.csv",df)
    return df
end

export MEME
function MEME(seqnames, seqs, treestring, outpath; verbosity=1, exports=true, code=MolecularEvolution.universal_code, optimize_branch_lengths=false, significance=0.05)
    exports && init_path(outpath)
    tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = MEME_global_fit(seqnames, seqs, treestring,
    verbosity=verbosity, code=code, optimize_branch_lengths=optimize_branch_lengths)
    alternative_params, alternative_LLs, null_params, null_LLs = MEME_hypotheses_fit(tree::FelNode, seqnames, seqs, GTRmat, F3x4_freqs, eq_freqs, code, verbosity = verbosity)
    LRTs, p_values = MEME_test(alternative_LLs, null_LLs, significance=significance)
    df = MEME_tabulate(LRTs, p_values, alternative_params, alternative_LLs, outpath, exports)
    #TODO: Branch test?
    return df
end
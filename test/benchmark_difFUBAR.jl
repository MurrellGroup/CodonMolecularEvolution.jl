function benchmark_grid_on_dataset(dir, versions_option, t)
    if any(endswith(".nex"), readdir(dir))
        nex_path = joinpath(dir, first(filter(endswith(".nex"), readdir(dir))))
        treestr, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(nex_path)
    elseif any(endswith(".nwk"), readdir(dir))
        nwk_path = joinpath(dir, first(filter(endswith(".nwk"), readdir(dir))))
        treestring_group_labeled, treestr, group_tags, tags, tag_colors = import_grouped_label_tree(nwk_path)
    end
    fasta_path = joinpath(dir, first(filter(endswith(".fasta"), readdir(dir))))
    analysis_name = ""

    seqnames, seqs = read_fasta(fasta_path)

    tree, tags, tag_colors, analysis_name = difFUBAR_init(analysis_name, treestr, tags, tag_colors, exports=false, verbosity=0)
    code = MolecularEvolution.universal_code
    tree, my_alpha, beta, GTRmat, F3x4_freqs, eq_freqs = difFUBAR_global_fit(seqnames, seqs, tree, generate_tag_stripper(tags), code, verbosity=0) #60s

    log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites = gridprep(tree, tags; verbosity = 1, foreground_grid = 6, background_grid = 4)
    heuristic_pick, nthreads = choose_grid_and_nthreads(tree, tags, num_groups, num_sites, alphagrid, omegagrid, background_omega_grid, code)
    if t > 0
        nthreads = t
    end
    @show typeof(heuristic_pick)
    grid_versions = [difFUBARBaseline(), difFUBARParallel(), difFUBARTreesurgery(), difFUBARTreesurgeryAndParallel()]
    option_map = Dict(1 => [heuristic_pick], 2 => [grid_versions[1]], 3 => [grid_versions[1], heuristic_pick], 4 => grid_versions)
    grid_versions_to_run = option_map[versions_option]
    num_taxa, num_sites, purity_ratio = length(getleaflist(tree)), tree.parent_message[1].partition.sites, get_purity_info(tree, tags, num_groups)[1]
    
    if dir == joinpath(@__DIR__, "data", "Ace2_no_background")
        println("Precompiling...")
        #Force compilation once
        for (i, grid) in enumerate(grid_versions_to_run)
            if grid isa difFUBARTreesurgery || grid isa difFUBARTreesurgeryAndParallel
                tree_local = deepcopy(tree)
            else
                tree_local = tree
            end
            log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites = gridprep(tree_local, tags; verbosity = 1, foreground_grid = 6, background_grid = 4)
            difFUBAR_grid(grid, tree_local, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites, nthreads; verbosity = 0, foreground_grid = 6, background_grid = 4)
        end
        GC.gc()
    end
    if versions_option == 4
        tree_tsp = deepcopy(tree)
    end
    con_lik_matrices = []
    times_and_bytes = []
    for (i, grid) in enumerate(grid_versions_to_run)
        if i == 4
            tree = tree_tsp
        end
        log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites = gridprep(tree, tags; verbosity = 1, foreground_grid = 6, background_grid = 4)
        (con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds), timed_grid, bytes_grid = @timed difFUBAR_grid(grid, tree, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites, nthreads; verbosity = 0, foreground_grid = 6, background_grid = 4)
        push!(con_lik_matrices, con_lik_matrix)
        push!(times_and_bytes, (timed_grid, bytes_grid))
    end
    if versions_option > 2
        for (con_lik_matrix, grid) in zip(con_lik_matrices[2:end], grid_versions_to_run[2:end])
            @assert isapprox(first(con_lik_matrices), con_lik_matrix) String(nameof(typeof(grid))) * " did not produce the same con_lik_matrix as the baseline version"
        end
    end

    f(time::Float64, bytes::Int64) = string(round(time, sigdigits=4)) * "s, " * string(round(bytes / 10^6, sigdigits=4)) * " M allocs"
    map(x -> f(x...), times_and_bytes), num_taxa, num_sites, purity_ratio, nthreads, typeof(heuristic_pick), grid_versions_to_run
end

function tabulate_global_fit_version(pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid)
    grid_size, num_sites = size(alloc_grid)

    r(s) = round(s,digits = 4);

    detected_sites = Int64[]
    detections = Vector{Float64}[] #legacy name - now includes all 4 "relevant" site posteriors
    param_means = Vector{Float64}[]

    ω1 = [c[2] for c in codon_param_vec];
    ω2 = [c[3] for c in codon_param_vec];
    alphas = [c[1] for c in codon_param_vec];
    ω1_greater_filt = ω1 .> ω2;
    ω2_greater_filt = ω2 .> ω1;
    ω1_pos_filt = ω1 .> 1.0;
    ω2_pos_filt = ω2 .> 1.0;

    for site in 1:num_sites
        ω1_greater_posterior = sum(alloc_grid[ω1_greater_filt,site])/sum(alloc_grid[:,site])
        ω2_greater_posterior = sum(alloc_grid[ω2_greater_filt,site])/sum(alloc_grid[:,site])
        ω1_pos_posterior = sum(alloc_grid[ω1_pos_filt,site])/sum(alloc_grid[:,site])
        ω2_pos_posterior = sum(alloc_grid[ω2_pos_filt,site])/sum(alloc_grid[:,site])
        detecs = [ω1_greater_posterior,ω2_greater_posterior,ω1_pos_posterior,ω2_pos_posterior]
        
        site_counts_ω1 = collapse_counts(ω1,alloc_grid[:,site], cases = omegagrid)
        site_counts_ω2 = collapse_counts(ω2,alloc_grid[:,site], cases = omegagrid)
        site_counts_alphas = collapse_counts(alphas,alloc_grid[:,site], cases = alphagrid)
        
        mean_alpha = sum(site_counts_alphas .* alphagrid)
        mean_ω1 = sum(site_counts_ω1 .* omegagrid)
        mean_ω2 = sum(site_counts_ω2 .* omegagrid)
        
        push!(detections,detecs)
        push!(param_means,[mean_alpha,mean_ω1,mean_ω2])
        
        if maximum(detecs) > pos_thresh
            push!(detected_sites,site)
        end
    end
    return detected_sites, detections, param_means, num_sites
end

function benchmark_global_fit_on_dataset(benchmark_name, dir, versions, nversions, exports, optimize_branch_lengths)
    if any(endswith(".nex"), readdir(dir))
        nex_path = joinpath(dir, first(filter(endswith(".nex"), readdir(dir))))
        treestr, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(nex_path)
    elseif any(endswith(".nwk"), readdir(dir))
        nwk_path = joinpath(dir, first(filter(endswith(".nwk"), readdir(dir))))
        treestring_group_labeled, treestr, group_tags, tags, tag_colors = import_grouped_label_tree(nwk_path)
    end
    fasta_path = joinpath(dir, first(filter(endswith(".fasta"), readdir(dir))))
    analysis_name = ""

    seqnames, seqs = read_fasta(fasta_path)

    tree, tags, tag_colors, analysis_name = difFUBAR_init(analysis_name, treestr, tags, tag_colors, exports=false, verbosity=0)
    code = MolecularEvolution.universal_code
    global_fits = versions
    trees = [tree, [deepcopy(tree) for _ = 1:(nversions - 1)]...]
    con_lik_matrices = []
    times_and_bytes = []
    if exports
        param_means_vec = []
        detections_vec = []
        detected_sites_vec = []
        pos_thresh = 0.95
        iters = 2500
        num_sites = 0
    end
    local_vars = []
    for (tree,global_fit) in zip(trees,global_fits)
        (tree, my_alpha, beta, GTRmat, F3x4_freqs, eq_freqs), timed_global_fit, bytes_global_fit = @timed global_fit(seqnames, seqs, tree, generate_tag_stripper(tags), code, verbosity=0, optimize_branch_lengths = optimize_branch_lengths) #60st

        push!(local_vars, (tree, GTRmat, F3x4_freqs))
        push!(times_and_bytes, (timed_global_fit, bytes_global_fit))
    end

    for (tree, GTRmat, F3x4_freqs) in local_vars
        log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites = gridprep(tree, tags; verbosity = 1, foreground_grid = 6, background_grid = 4)
        heuristic_pick, nthreads = choose_grid_and_nthreads(tree, tags, num_groups, num_sites, alphagrid, omegagrid, background_omega_grid, code)
        con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = difFUBAR_grid(heuristic_pick, tree, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites, nthreads; verbosity = 0, foreground_grid = 6, background_grid = 4)
        if exports
            alloc_grid,theta = difFUBAR_sample(con_lik_matrix, iters, verbosity = false)
            detected_sites, detections, param_means, ns = tabulate_global_fit_version(pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid)
            num_sites = ns
            push!(param_means_vec, param_means)
            push!(detections_vec, detections)
            push!(detected_sites_vec, detected_sites)
        end
        push!(con_lik_matrices, con_lik_matrix)
    end
    if exports
        plots = []
        titles = String[]
        #Plot means
        categories = ["α", "ω1", "ω2"]
        for (i, category) in enumerate(categories)
            mean_codon = [p[i] for p in param_means_vec[1]]
            mean_nuc_codon = [p[i] for p in param_means_vec[2]]
            x = 0.0:0.01:maximum(mean_codon)
            p = plot(x, x, color=:black, label=false, aspect_ratio=:equal)
            scatter!( mean_codon, mean_nuc_codon, color=:blue, label=false, aspect_ratio=:equal)
            push!(titles, category * " means")
            push!(plots, p)
        end

        #Plot posterior probabilities
        categories = ["P(ω1 > ω2)", "P(ω2 > ω1)", "P(ω1 > 1)", "P(ω2 > 1)"]
        p2 = plot(layout = (1, length(categories)), thickness_scaling = 0.5)
        for (i, category) in enumerate(categories)
            codon = [p[i] for p in detections_vec[1]]
            nuc_codon = [p[i] for p in detections_vec[2]]
            x = 0.0:0.01:1.0
            p = plot(x, x, color=:black, label=false, aspect_ratio=:equal)
            scatter!(codon, nuc_codon, color=:red, label=false, aspect_ratio=:equal)
            push!(titles, category * " posterior")
            push!(plots, p)
        end
        l = @layout([° ° ° _; ° ° ° °])
        plt = plot(plots..., layout = l, size = (800, 400), title=reshape(titles, 1, 7), plot_title="Codon vs. Nucleotide+Codon model fit", xlabel="Codon", ylabel="Nucleotide+Codon", left_margin=50px, bottom_margin=50px, thickness_scaling = 0.5)
        savefig(benchmark_name * "_means_and_posteriors.pdf")
    end

    f(time::Float64, bytes::Int64) = string(round(time, sigdigits=4)) * "s, " * string(round(bytes / 10^6, sigdigits=4)) * " M allocs"
    return map(x -> f(x...), times_and_bytes), map(x -> maximum(abs.(con_lik_matrices[1] - x)), con_lik_matrices[2:end]), map(x -> sum(abs.(con_lik_matrices[1] - x)) / length(con_lik_matrices[1]), con_lik_matrices[2:end])
end
"""
    CodonMolecularEvolution.benchmark_grid(benchmark_name; exports=true, versions_option=1, t::Integer=0, data=1:5)
Benchmarks different implementations of the difFUBAR_grid algorithm. Results of the benchmark are printed out as a DataFrame and saved to a CSV file.
- `benchmark_name` is the filepath to where the benchmark will be saved, if exports
- `versions_option` have 4 different options:
    - 1. default option, only run heuristic top pick
    - 2. only run baseline version
    - 3. run heuristic top pick and baseline version
    - 4. run all versions
- `t` is the number of threads you want to use in the parallel versions. If specified and non-zero, this will override the number of threads chosen by the heuristic.
- `data` is the range/vector of datasets to run the benchmark on. By default, this is 1:5. These are the enumerated datasets:
    - 1. Ace2nobackground
    - 2. Ace2reallytiny
    - 3. Ace2tiny
    - 4. ParvoVP
    - 5. ParvoVPregrouped
"""
function benchmark_grid(benchmark_name; exports=true, versions_option=1, t::Integer=0, data=1:5)
    #Create the export directory, if required
    splt = splitpath(benchmark_name)[1:end-1]
    if length(splt) > 0
        exports && mkpath(joinpath(splt))
    end
    nversions = 4
    data_dir = joinpath(@__DIR__, "data")
    n = length(readdir(data_dir))
    timings = fill("", n, nversions)
    num_sites_vec = Vector{Int64}(undef, n)
    num_taxa_vec = Vector{Int64}(undef, n)
    datasets = Vector{String}(undef, n)
    purity_ratio_vec = Vector{Float64}(undef, n)
    nthreads_vec = Vector{Int64}(undef, n)
    heuristic_picks = Vector{DataType}(undef, n)

    version_map = Dict(difFUBARBaseline => 1, difFUBARParallel => 2, difFUBARTreesurgery => 3, difFUBARTreesurgeryAndParallel => 4)

    for (i, dataset) in enumerate(readdir(data_dir))
        if !(i in data)
            continue
        end
        @show dataset
        dir = joinpath(data_dir, dataset)
        timing, num_taxa, num_sites, purity_ratio, nthreads, heuristic_pick, versions_ran = benchmark_grid_on_dataset(dir, versions_option, t)
        for (time_and_bytes, version) in zip(timing, versions_ran)
            timings[i, version_map[typeof(version)]] = time_and_bytes
        end

        num_taxa_vec[i] = num_taxa
        num_sites_vec[i] = num_sites
        datasets[i] = dataset
        purity_ratio_vec[i] = purity_ratio
        nthreads_vec[i] = nthreads
        heuristic_picks[i] = heuristic_pick
    end

    df = DataFrame(hcat(datasets, num_taxa_vec, num_sites_vec, purity_ratio_vec, nthreads_vec, heuristic_picks, timings)[data, :], [:dataset, :num_taxa, :num_sites, :purity_ratio, :nthreads, :heuristic_pick, :baseline, :parallel, :treesurgery, :treesurgery_and_parallel])
    println(df)
    exports && CSV.write(benchmark_name*".csv", df);
end

"""
    CodonMolecularEvolution.benchmark_global_fit(benchmark_name; exports=true, data=1:5, optimize_branch_lengths=false)
Benchmarks different implementations of the difFUBAR_global_fit algorithm. Results of the benchmark are printed out as a DataFrame and saved to a CSV file. 
Uses the heuristic top pick to generate con lik matrices. Compares difference in con lik matrices.
If exports is true, it also runs MCMCs on the con lik matrices and plots the means and posteriors of the different versions.
- `benchmark_name` is the filepath to where the benchmark will be saved, if exports.
- `data` is the range/vector of datasets to run the benchmark on. By default, this is 1:5. These are the enumerated datasets:
    - 1. Ace2nobackground
    - 2. Ace2reallytiny
    - 3. Ace2tiny
    - 4. ParvoVP
    - 5. ParvoVPregrouped
- `optimize_branch_lengths` is an option that can be either `true`, `false` or `"detect"`
"""
function benchmark_global_fit(benchmark_name; exports=true, data=1:5, optimize_branch_lengths=false)
    #Create the export directory, if required
    splt = splitpath(benchmark_name)[1:end-1]
    if length(splt) > 0
        exports && mkpath(joinpath(splt))
    end
    versions = [difFUBAR_global_fit, difFUBAR_global_fit_2steps]
    nversions = length(versions)
    data_dir = joinpath(@__DIR__, "data")
    n = length(readdir(data_dir))
    timings = fill("", n, nversions)
    max_diffs = Vector{Vector{Float64}}(undef, n)
    avg_diffs = Vector{Vector{Float64}}(undef, n)
    datasets = Vector{String}(undef, n)

    for (i, dataset) in enumerate(readdir(data_dir))
        if !(i in data)
            continue
        end
        @show dataset
        dir = joinpath(data_dir, dataset)
        timing, max_diff, avg_diff = benchmark_global_fit_on_dataset(benchmark_name*dataset, dir, versions, nversions, exports, optimize_branch_lengths)
        timings[i, :] .= timing
        datasets[i] = dataset
        max_diffs[i] = max_diff
        avg_diffs[i] = avg_diff
    end

    df = DataFrame(hcat(datasets, timings, max_diffs, avg_diffs)[data, :], [:dataset, :global_fit, :global_fit_2steps, :max_diff, :avg_diff])
    println(df)
    exports && CSV.write(benchmark_name*".csv", df);
end

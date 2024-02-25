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
    @show nameof(heuristic_pick)
    grid_versions = [difFUBAR_grid_baseline, difFUBAR_grid_parallel, difFUBAR_grid_treesurgery, difFUBAR_grid_treesurgery_and_parallel]
    option_map = Dict(1 => [heuristic_pick], 2 => [grid_versions[1]], 3 => [grid_versions[1], heuristic_pick], 4 => grid_versions)
    grid_versions_to_run = option_map[versions_option]
    num_taxa, num_sites, purity_ratio = length(getleaflist(tree)), tree.message[1].sites, get_purity_info(tree, tags, num_groups)[1]
    
    if dir == joinpath(@__DIR__, "data", "Ace2_no_background")
        println("Precompiling...")
        #Force compilation once
        for (i, grid) in enumerate(grid_versions_to_run)
            if grid == difFUBAR_grid_treesurgery || grid == difFUBAR_grid_treesurgery_and_parallel
                tree_local = deepcopy(tree)
            else
                tree_local = tree
            end
            log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites = gridprep(tree_local, tags; verbosity = 1, foreground_grid = 6, background_grid = 4)
            grid(tree_local, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites, nthreads; verbosity = 0, foreground_grid = 6, background_grid = 4)
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
        (con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds), timed_grid, bytes_grid = @timed grid(tree, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites, nthreads; verbosity = 0, foreground_grid = 6, background_grid = 4)
        push!(con_lik_matrices, con_lik_matrix)
        push!(times_and_bytes, (timed_grid, bytes_grid))
    end
    if versions_option > 2
        for (con_lik_matrix, grid) in zip(con_lik_matrices[2:end], grid_versions_to_run[2:end])
            @assert isapprox(first(con_lik_matrices), con_lik_matrix) String(nameof(grid)) * " did not produce the same con_lik_matrix as the baseline version"
        end
    end

    f(time::Float64, bytes::Int64) = string(round(time, sigdigits=4)) * "s, " * string(round(bytes / 10^6, sigdigits=4)) * " M allocs"
    map(x -> f(x...), times_and_bytes), num_taxa, num_sites, purity_ratio, nthreads, nameof(heuristic_pick), grid_versions_to_run
end

function benchmark_global_fit_on_dataset(dir)
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
    global_fits = [difFUBAR_global_fit, difFUBAR_global_fit_2steps]
    trees = [tree, deepcopy(tree)]
    con_lik_matrices = []
    times_and_bytes = []
    local_vars = []
    for (tree,global_fit) in zip(trees,global_fits)
        (tree, my_alpha, beta, GTRmat, F3x4_freqs, eq_freqs), timed_global_fit, bytes_global_fit = @timed global_fit(seqnames, seqs, tree, generate_tag_stripper(tags), code, verbosity=0) #60s
        @show GTRmat

        push!(local_vars, (tree, GTRmat, F3x4_freqs))
        push!(times_and_bytes, (timed_global_fit, bytes_global_fit))
    end
    
    for (tree, GTRmat, F3x4_freqs) in local_vars
        log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites = gridprep(tree, tags; verbosity = 1, foreground_grid = 6, background_grid = 4)
        heuristic_pick, nthreads = choose_grid_and_nthreads(tree, tags, num_groups, num_sites, alphagrid, omegagrid, background_omega_grid, code)
        con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = heuristic_pick(tree, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites, nthreads; verbosity = 0, foreground_grid = 6, background_grid = 4)
        push!(con_lik_matrices, con_lik_matrix)
    end

    f(time::Float64, bytes::Int64) = string(round(time, sigdigits=4)) * "s, " * string(round(bytes / 10^6, sigdigits=4)) * " M allocs"
    return map(x -> f(x...), times_and_bytes), maximum(abs.(con_lik_matrices[1] - con_lik_matrices[2])), sum(abs.(con_lik_matrices[1] - con_lik_matrices[2])) / length(con_lik_matrices[1])
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
    heuristic_picks = Vector{Symbol}(undef, n)

    version_map = Dict(difFUBAR_grid_baseline => 1, difFUBAR_grid_parallel => 2, difFUBAR_grid_treesurgery => 3, difFUBAR_grid_treesurgery_and_parallel => 4)

    for (i, dataset) in enumerate(readdir(data_dir))
        if !(i in data)
            continue
        end
        @show dataset
        dir = joinpath(data_dir, dataset)
        timing, num_taxa, num_sites, purity_ratio, nthreads, heuristic_pick, versions_ran = benchmark_grid_on_dataset(dir, versions_option, t)
        for (time_and_bytes, version) in zip(timing, versions_ran)
            timings[i, version_map[version]] = time_and_bytes
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
    CodonMolecularEvolution.benchmark_global_fit(benchmark_name; exports=true, data=1:5)
Benchmarks different implementations of the difFUBAR_global_fit algorithm. Results of the benchmark are printed out as a DataFrame and saved to a CSV file. 
Uses the heuristic top pick to generate con lik matrices. Compares difference in con lik matrices.
- `benchmark_name` is the filepath to where the benchmark will be saved, if exports
- `data` is the range/vector of datasets to run the benchmark on. By default, this is 1:5. These are the enumerated datasets:
    - 1. Ace2nobackground
    - 2. Ace2reallytiny
    - 3. Ace2tiny
    - 4. ParvoVP
    - 5. ParvoVPregrouped
"""
function benchmark_global_fit(benchmark_name; exports=true, data=1:5)
    #Create the export directory, if required
    splt = splitpath(benchmark_name)[1:end-1]
    if length(splt) > 0
        exports && mkpath(joinpath(splt))
    end
    nversions = 2
    data_dir = joinpath(@__DIR__, "data")
    n = length(readdir(data_dir))
    timings = fill("", n, nversions)
    max_diffs = Vector{Float64}(undef, n)
    avg_diffs = Vector{Float64}(undef, n)
    datasets = Vector{String}(undef, n)

    version_map = Dict(difFUBAR_global_fit => 1, difFUBAR_global_fit_2steps => 2)

    for (i, dataset) in enumerate(readdir(data_dir))
        if !(i in data)
            continue
        end
        @show dataset
        dir = joinpath(data_dir, dataset)
        timing, max_diff, avg_diff = benchmark_global_fit_on_dataset(dir)
        timings[i, :] .= timing
        datasets[i] = dataset
        max_diffs[i] = max_diff
        avg_diffs[i] = avg_diff
    end

    df = DataFrame(hcat(datasets, timings, max_diffs, avg_diffs)[data, :], [:dataset, :global_fit, :global_fit_2steps, :max_diff, :avg_diff])
    println(df)
    exports && CSV.write(benchmark_name*".csv", df);
end

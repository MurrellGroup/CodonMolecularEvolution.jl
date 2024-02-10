using MolecularEvolution, CodonMolecularEvolution, Glob, BenchmarkTools, CSV, DataFrames, LinearAlgebra

BLAS.set_num_threads(1)

function test_grids(dir)
    @show dir
    nex_path = glob(dir*"*.nex")[1]
    fasta_path = glob(dir*"*.fasta")[1]
    analysis_name = ""

    treestr, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(nex_path)
    seqnames, seqs = read_fasta(fasta_path)

    tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestr, tags, tag_colors, exports=false, verbosity=1)
    code = MolecularEvolution.universal_code
    tree, my_alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #60s

    if dir == "./test/data/Ace2_no_background/"
        println("Precompiling...")
        tree_ts = deepcopy(tree)
        tree_tsp = deepcopy(tree)

        #Force compilation once
        log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites = CodonMolecularEvolution.gridprep(tree, tags; verbosity = 1, foreground_grid = 6, background_grid = 4)
        CodonMolecularEvolution.difFUBAR_grid_parallel(tree, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites; verbosity = 0, foreground_grid = 6, background_grid = 4)
    
        log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites = CodonMolecularEvolution.gridprep(tree_ts, tags; verbosity = 1, foreground_grid = 6, background_grid = 4)
        CodonMolecularEvolution.difFUBAR_grid_treesurgery(tree_ts, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites; verbosity = 0, foreground_grid = 6, background_grid = 4)

        log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites = CodonMolecularEvolution.gridprep(tree_tsp, tags; verbosity = 1, foreground_grid = 6, background_grid = 4)
        CodonMolecularEvolution.difFUBAR_grid_treesurgery_and_parallel(tree_tsp, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites; verbosity = 0, foreground_grid = 6, background_grid = 4)
        GC.gc()
    end
    tree_ts = deepcopy(tree)
    tree_tsp = deepcopy(tree)

    #Benchmark
    #Warning: the benchmarks for the tree_surgery methods are not accurate because the tree gets truncated after first use
    #Consider using @timed instead
    (con_lik_matrix_no_prune, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds), timed_difFUBAR, bytes_difFUBAR = @timed CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=0, foreground_grid=6, background_grid=4) # 50s
    
    log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites = CodonMolecularEvolution.gridprep(tree, tags; verbosity = 1, foreground_grid = 6, background_grid = 4)
    (con_lik_matrix_parallel, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds), timed_difFUBAR_parallel, bytes_difFUBAR_parallel = @timed CodonMolecularEvolution.difFUBAR_grid_parallel(tree, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites; verbosity = 0, foreground_grid = 6, background_grid = 4)
    
    log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites = CodonMolecularEvolution.gridprep(tree_ts, tags; verbosity = 1, foreground_grid = 6, background_grid = 4)
    (con_lik_matrix_treesurgery, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds), timed_difFUBAR_treesurgery, bytes_difFUBAR_treesurgery = @timed CodonMolecularEvolution.difFUBAR_grid_treesurgery(tree_ts, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites; verbosity = 0, foreground_grid = 6, background_grid = 4)

    log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites = CodonMolecularEvolution.gridprep(tree_tsp, tags; verbosity = 1, foreground_grid = 6, background_grid = 4)
    (con_lik_matrix_treesurgery_and_parallel, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds), timed_difFUBAR_treesurgery_and_parallel, bytes_difFUBAR_treesurgery_and_parallel = @timed CodonMolecularEvolution.difFUBAR_grid_treesurgery_and_parallel(tree_tsp, tags, GTRmat, F3x4_freqs, code, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, background_omega_grid, param_kinds, is_background, num_groups, num_sites; verbosity = 0, foreground_grid = 6, background_grid = 4)
    
    @assert con_lik_matrix_no_prune == con_lik_matrix_parallel
    @assert con_lik_matrix_no_prune == con_lik_matrix_treesurgery
    @assert con_lik_matrix_no_prune == con_lik_matrix_treesurgery_and_parallel

    f(time::Float64, bytes::Int64) = string(round(time, digits=4)) * "s, " * string(round(bytes / 10^6, digits=4)) * " M allocs"
    [f(timed_difFUBAR, bytes_difFUBAR), f(timed_difFUBAR_parallel, bytes_difFUBAR_parallel), f(timed_difFUBAR_treesurgery, bytes_difFUBAR_treesurgery), f(timed_difFUBAR_treesurgery_and_parallel, bytes_difFUBAR_treesurgery_and_parallel)], length(getleaflist(tree)), tree.message[1].sites
end

dir_outer = "./test/data"
n = length(readdir(dir_outer))
timings = Matrix{String}(undef, n, 4)
num_sites_vec = Vector{Int64}(undef, n)
num_taxa_vec = Vector{Int64}(undef, n)
dir_names = Vector{String}(undef, n)

for (i, dir_name) in enumerate(readdir(dir_outer))
    dir = joinpath(dir_outer, dir_name, "")
    timing, num_taxa, num_sites = test_grids(dir)
    timings[i, :] = timing
    num_taxa_vec[i] = num_taxa
    num_sites_vec[i] = num_sites
    dir_names[i] = dir_name
end

df = DataFrame(hcat(dir_names, num_taxa_vec, num_sites_vec, timings), [:test_case, :num_taxa, :num_sites, :difFUBAR, :difFUBAR_parallel, :difFUBAR_treesurgery, :difFUBAR_treesurgery_and_parallel])
CSV.write("./test/timings-t16-BLAS-1-ParvoVP-included.csv", df) 
exit()

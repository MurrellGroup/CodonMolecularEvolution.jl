using Revise

using MolecularEvolution, CodonMolecularEvolution

treestr, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(".\\Trees\\Ace2_tiny_tagged.nex")
seqnames, seqs = read_fasta(".\\Trees\\Ace2_tiny_tagged.fasta")
analysis_name = ".\\Paralllel\\Results\\"

tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestr, tags, tag_colors, exports=false, verbosity=1)
code = MolecularEvolution.universal_code
tree, my_alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #60s

tree_no_prune = deepcopy(tree)
tree_prune_1 = deepcopy(tree)
tree_prune_2 = deepcopy(tree)

CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=0, foreground_grid=6, background_grid=4)
time_difFUBAR = @elapsed con_lik_matrix_no_prune, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid(tree_no_prune, tags, GTRmat, F3x4_freqs, code, verbosity=0, foreground_grid=6, background_grid=4) # 50s
time_difFUBAR_parallel = @elapsed con_lik_matrix_parallel, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_parallel(tree_prune_1, tags, GTRmat, F3x4_freqs, code, verbosity=0, foreground_grid=6, background_grid=4) # 28s
@assert con_lik_matrix_no_prune == con_lik_matrix_parallel
#time_difFUBAR_prune_max = @elapsed con_lik_matrix_prune_2, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_1(tree_prune_2, tags, GTRmat, F3x4_freqs, code, verbosity=0, foreground_grid=6, background_grid=4) # 28s
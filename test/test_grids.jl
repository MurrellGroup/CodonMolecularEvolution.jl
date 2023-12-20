using Revise

using MolecularEvolution, CodonMolecularEvolution

# On Windows
nex_path = ".\\Trees\\Ace2_reallytiny_tagged.nex"
fasta_path = ".\\Trees\\Ace2_reallytiny_tagged.fasta"
analysis_name = ".\\Paralllel\\Results\\"

# On macOS
nex_path = "./TreeSurgery/Trees/RealData/Ace2_reallytiny_tagged.nex"
fasta_path = "./TreeSurgery/Trees/RealData/Ace2_reallytiny_tagged.fasta"
analysis_name = "./TreeSurgery/Results/"

treestr, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(nex_path)
seqnames, seqs = read_fasta(fasta_path)

tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestr, tags, tag_colors, exports=false, verbosity=1)
code = MolecularEvolution.universal_code
tree, my_alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #60s

tree_no_prune = deepcopy(tree)
tree_parallel = deepcopy(tree)
tree_chunks = deepcopy(tree)
tree_maxi = deepcopy(tree)

CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=0, foreground_grid=6, background_grid=4)
time_difFUBAR = @elapsed con_lik_matrix_no_prune, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid(tree_no_prune, tags, GTRmat, F3x4_freqs, code, verbosity=0, foreground_grid=6, background_grid=4) # 50s
time_difFUBAR_parallel = @elapsed con_lik_matrix_parallel, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_parallel(tree_parallel, tags, GTRmat, F3x4_freqs, code, verbosity=0, foreground_grid=6, background_grid=4) # 28s
time_difFUBAR_chunks = @elapsed con_lik_matrix_chunks, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_chunks(tree_chunks, tags, GTRmat, F3x4_freqs, code, verbosity=0, foreground_grid=6, background_grid=4)
time_difFUBAR_maxi = @elapsed con_lik_matrix_maxi, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_maxi(tree_maxi, tags, GTRmat, F3x4_freqs, code, verbosity=0, foreground_grid=6, background_grid=4) # 28s
@assert con_lik_matrix_no_prune == con_lik_matrix_parallel
@assert con_lik_matrix_no_prune == con_lik_matrix_chunks
@assert con_lik_matrix_no_prune == con_lik_matrix_maxi

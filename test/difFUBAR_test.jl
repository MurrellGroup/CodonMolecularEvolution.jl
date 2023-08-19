analysis_name = "nobackground/Ace2"
seqnames,seqs = read_fasta("data/Ace2_tiny_test.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/Ace2_no_background.nex")
df,results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports = false, verbosity = 0)
@assert size(df) == (19, 8)
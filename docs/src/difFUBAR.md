# difFUBAR
An interface to [difFUBAR (awaiting the correct link)](https://academic.oup.com/mbe/article/30/5/1196/998247)
## Quick example
Reads codon sequences from [this FASTA file (awaiting the correct link)](https://raw.githubusercontent.com/MurrellGroup/MolecularEvolution.jl/main/docs/src/Flu.fasta), and a phylogeny from [this Newick tree file (awaiting the correct link)](https://raw.githubusercontent.com/MurrellGroup/MolecularEvolution.jl/main/docs/src/Flu.tre).
```julia
using MolecularEvolution, FASTX, CodonMolecularEvolution
analysis_name = "output/Ace2"
seqnames,seqs = read_fasta("data/Ace2_tiny_test.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/Ace2_no_background.nex")
df,results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name);
```
```
Step 1: Initialization. If exports = true, tree showing the assignment of branches to groups/colors will be exported to: output/Ace2_tagged_input_tree.svg.
Step 2: Optimizing global codon model parameters.
0.0% 29.0% 58.0% 87.0% 
Step 4: Running Gibbs sampler to infer site categories.
Step 5: Tabulating and plotting. Detected sites:
Site 3 - P(ω1 > ω2):0.006; P(ω2 > ω1):0.981; P(ω1 > 1):0.113; P(ω2 > 1):0.9755
Site 13 - P(ω1 > ω2):0.0205; P(ω2 > ω1):0.9565; P(ω1 > 1):0.14; P(ω2 > 1):0.9175

If exports = true, writing results for all sites to CSV: output/Ace2_posteriors.csv
Plotting alpha and omega distributions. If exports = true, saved as output/Ace2_violin_*.pdf
```
TODO: add output figures

## Interface
```@docs
difFUBAR
```

### For a simple and often optimal configuration
- Launch Julia in the following manner: `julia -t auto`
- Keep the default values of the kwargs `version` and `t`
This lets Julia decide the amount of Julia threads and lets CodonMolecularEvolution.jl decide the [difFUBARGrid](@ref) subtype to dispatch on and the degree of parallelization.

## difFUBARGrid
**Subtypes** that decide which method to use for the grid likelihood computations.
```@docs
difFUBARBaseline
difFUBARParallel
difFUBARTreesurgery
difFUBARTreesurgeryAndParallel
```
<img src="https://user-images.githubusercontent.com/1152087/188331266-5e03565b-00a7-490c-a616-50598ca46010.png" width="140">

# CodonMolecularEvolution.jl

<!---[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MurrellGroup.github.io/CodonMolecularEvolution.jl/stable/)--->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MurrellGroup.github.io/CodonMolecularEvolution.jl/dev/)
[![Coverage](https://codecov.io/gh/MurrellGroup/CodonMolecularEvolution.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/CodonMolecularEvolution.jl)
[![Build Status](https://github.com/MurrellGroup/CodonMolecularEvolution.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/CodonMolecularEvolution.jl/actions/workflows/CI.yml?query=branch%3Amain)

Descendant of [MolecularEvolution.jl](https://github.com/MurrellGroup/MolecularEvolution.jl), specializing in codon models. The current dN/dS catalogue consists of

- FUBAR
- difFUBAR

To run our models in the cloud, see [ColabMolecularEvolution](https://github.com/MurrellGroup/ColabMolecularEvolution).
## Installation
In a julia REPL, type
`]add CodonMolecularEvolution`

## Documentation
https://MurrellGroup.github.io/CodonMolecularEvolution.jl

## Extensions
To trigger our plotting extension, one must `using Plots, Phylo`.
## Examples
This section demonstrates how to use our dN/dS models with real biological data, showcasing both the input requirements and the type of evolutionary insights you can obtain.
```julia
using MolecularEvolution, FASTX, CodonMolecularEvolution
```
**Note**: To enable plotting, one would add `using Plots, Phylo`.
### FUBAR
<a href="https://academic.oup.com/mbe/article/30/5/1196/998247">
    <img src="https://img.shields.io/badge/DOI-10.1093%2Fmolbev%2Fmst030-blue.svg" alt="DOI"/>
</a>

We perform FUBAR analysis on [this FASTA file](https://raw.githubusercontent.com/MurrellGroup/CodonMolecularEvolution.jl/main/test/data/Ace2_with_bat/Ace2_with_bat.fasta), and a phylogeny from [this Newick tree file](https://raw.githubusercontent.com/MurrellGroup/CodonMolecularEvolution.jl/main/test/data/Ace2_with_bat/Ace2_with_bat.tre).
```julia
# Read data
outdir = "fubar"
seqnames, seqs = read_fasta("Ace2_with_bat.fasta")
treestring = readlines("Ace2_with_bat.tre")[1]
# Perform analysis
fgrid = alphabetagrid(seqnames, seqs, treestring);
method = DirichletFUBAR() # Standard FUBAR
fubar_df, params = FUBAR_analysis(method, fgrid, analysis_name=outdir*"/Ace2_with_bat");
```
The above call will write results to `outdir*"/Ace2_with_bat"`, but here we only display what's written to sdout
```
Step 1: Initialization.
Step 2: Optimizing global codon model parameters.
Optimized single α,β LL=-66072.56532041529 with α=1.940881909909932 and β=0.7757100788818039.
Step 3: Calculating conditional likelihoods.
Step 4: Model fitting.
Step 5: Tabulating results.
19 sites with positive selection above threshold.
487 sites with purifying selection above threshold.
```
Now we show some of the above-threshold sites
```julia
above_threshold_sites = fubar_df[(fubar_df.positive_posterior .> 0.95) .| (fubar_df.purifying_posterior .> 0.95), :]
above_threshold_sites[1:10, :]
```
```
10×5 DataFrame
 Row │ site   positive_posterior  purifying_posterior  beta_posterior_mean  alpha_pos_mean 
     │ Int64  Float64             Float64              Float64              Float64        
─────┼─────────────────────────────────────────────────────────────────────────────────────
   1 │     1         0.00564089            0.985948              0.0126755        1.42783
   2 │     4         5.12413e-5            0.998737              0.218888         1.12801
   3 │     6         0.0092943             0.981231              0.0227715        1.47458
   4 │     9         0.0118501             0.964818              0.0658507        0.392239
   5 │    11         2.29686e-8            0.999994              0.277966         2.10407
   6 │    15         1.43323e-6            0.999887              0.20953          1.42799
   7 │    17         0.000436892           0.995299              0.199044         0.917965
   8 │    19         8.04229e-9            0.999999              0.152897         2.47035
   9 │    25         0.000141827           0.996911              0.275294         1.36792
  10 │    27         0.995895              3.74494e-5            2.94698          1.05116
```
### difFUBAR (awaiting the correct link)
<a href="https://academic.oup.com/mbe/article/30/5/1196/998247">
    <img src="https://img.shields.io/badge/DOI-10.1093%2Fmolbev%2Fmst030-blue.svg" alt="DOI"/>
</a>

We perform difFUBAR analysis on [this FASTA file](https://raw.githubusercontent.com/MurrellGroup/CodonMolecularEvolution.jl/main/test/data/Ace2_no_background/Ace2_tiny_test.fasta), and a tagged phylogeny from [this NEXUS tree file](https://raw.githubusercontent.com/MurrellGroup/CodonMolecularEvolution.jl/main/test/data/Ace2_no_background/Ace2_no_background.nex).
```julia
# Read data
analysis_name = "difFUBAR/Ace2_tiny"
seqnames,seqs = read_fasta("data/Ace2_tiny_test.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/Ace2_no_background.nex")
# Perform analysis
df,results = difFUBAR(seqnames, seqs, treestring, tags, analysis_name);
```
`difFUBAR` will write its results to the export directory `analysis_name`,
but here we show a short-hand representation of the results, which are written to stdout.

```
Step 1: Initialization. If exports = true, tree showing the assignment of branches to groups/colors will be exported to: output/Ace2_tagged_input_tree.svg.
Step 2: Optimizing global codon model parameters.
0.0% 29.0% 58.0% 87.0% 
Step 4: Running Gibbs sampler to infer site categories.
Step 5: Tabulating and plotting. Detected sites:
Site 3 - P(ω1 > ω2):0.0075; P(ω2 > ω1):0.9805; P(ω1 > 1):0.1205; P(ω2 > 1):0.9675
Site 13 - P(ω1 > ω2):0.0175; P(ω2 > ω1):0.9605; P(ω1 > 1):0.13; P(ω2 > 1):0.9305
```

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
## Examples
Throughout, we will read codon sequences from [this FASTA file](https://raw.githubusercontent.com/MurrellGroup/CodonMolecularEvolution.jl/main/test/data/Ace2_no_background/Ace2_tiny_test.fasta), and a tagged phylogeny from [this NEXUS tree file](https://raw.githubusercontent.com/MurrellGroup/CodonMolecularEvolution.jl/main/test/data/Ace2_no_background/Ace2_no_background.nex)
```julia
using MolecularEvolution, FASTX, CodonMolecularEvolution
analysis_name = "output/Ace2"
seqnames,seqs = read_fasta("data/Ace2_tiny_test.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/Ace2_no_background.nex")
```
### FUBAR
<a href="https://academic.oup.com/mbe/article/30/5/1196/998247">
    <img src="https://img.shields.io/badge/DOI-10.1093%2Fmolbev%2Fmst030-blue.svg" alt="DOI"/>
</a>

We perform FUBAR analysis on the data.
```julia
# In FUBAR, we treat the whole alignment homogeneously and don't care about the tags
stripped_treestring = CodonMolecularEvolution.generate_tag_stripper(tags)(treestring)
fgrid = alphabetagrid(seqnames, seqs, stripped_treestring);
fubar_df, params = FUBAR(fgrid, outdir*"/fubarACE2");
```
The above call will write results to `outdir*"/fubarACE2"`, but here we only display what's written to sdout
```
Step 1: Initialization.
Step 2: Optimizing global codon model parameters.
Optimized single α,β LL=-357.4118275678246 with α=1.2245983025103175 and β=0.9236464501066286.
Step 3: Calculating conditional likelihoods.
Step 4: Model fitting.
Step 5: Tabulating results and saving plots.
0 sites with positive selection above threshold.
2 sites with purifying selection above threshold.
```
Now we show the above-threshold sites:
```julia
significant_sites = fubar_df[(fubar_df.positive_posterior .> 0.95) .| (fubar_df.purifying_posterior .> 0.95), :]
```
```
2×5 DataFrame
 Row │ site   positive_posterior  purifying_posterior  beta_pos_mean  alpha_pos_mean 
     │ Int64  Float64             Float64              Float64        Float64        
─────┼───────────────────────────────────────────────────────────────────────────────
   1 │    11           0.0185862             0.953521       0.852326         4.26415
   2 │    16           0.0171059             0.962277       0.48187          3.4233
```
### difFUBAR (awaiting the correct link)
<a href="https://academic.oup.com/mbe/article/30/5/1196/998247">
    <img src="https://img.shields.io/badge/DOI-10.1093%2Fmolbev%2Fmst030-blue.svg" alt="DOI"/>
</a>

We perform difFUBAR analysis on the data.
```julia
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
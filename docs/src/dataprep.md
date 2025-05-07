# Data preparation

## Branch tagging

Starting from this phylogeny [untagged newick file](https://raw.githubusercontent.com/MurrellGroup/CodonMolecularEvolution.jl/main/docs/src/untagged_example.nwk), 
![](figures/untagged_example_tagged_input_tree.svg)

an alternative we can use for branch tagging is the [phylotree web application](https://phylotree.hyphy.org/). See [phylotree tutorial](https://hyphy.org/tutorials/phylotree/) for details on how to do this.

Then, after a few button clicks in the app, we get a [tagged newick file](https://raw.githubusercontent.com/MurrellGroup/CodonMolecularEvolution.jl/main/docs/src/tagged_example.nwk)

![](figures/tagged_example_tagged_input_tree.svg)

We can extract the `treestring` and the `tags` from the newick file with [`import_labeled_phylotree_newick`](@ref)
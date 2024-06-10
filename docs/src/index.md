```@meta
CurrentModule = CodonMolecularEvolution
```

# CodonMolecularEvolution

Documentation for [CodonMolecularEvolution](https://github.com/MurrellGroup/CodonMolecularEvolution.jl).

Descendant of [MolecularEvolution.jl](https://github.com/MurrellGroup/MolecularEvolution.jl), specializing in codon models.

### A Julia package that contains a collection of popular/(state of the art) codon model methods.

CodonMolecularEvolution.jl leverages the MolecularEvolution.jl framework to provide a performant interface to some codon model methods that are used for phylogenetic analysis.

### Collection of codon model methods
- [difFUBAR](@ref): Scalable Bayesian comparison of selection pressure
    - Perform a site-wise comparison of evolutionary pressure between two selected sets of branches.

### Design principles
- Specificity
    - We keep our implementations on a high level of abstraction and let `MolecularEvolution.jl` take care of the low-level behaviour.
- Scalability
    - Our exported codon model methods should scale to large, real-world datasets. To keep the memory footprint down, we use `MolecularEvolution.jl`s `LazyPartition` when possible.
- Performance
    - While scalability takes precedence over speed?, we try to maintain competitive runtimes by using e.g. computational shortcuts and parallelization whenever we can.

### Authors
Ben decides what to put here.

```@index
```

```@autodocs
Modules = [CodonMolecularEvolution]
```

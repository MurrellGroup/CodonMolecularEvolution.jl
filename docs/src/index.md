```@meta
CurrentModule = CodonMolecularEvolution
```

# CodonMolecularEvolution

### A Julia package for popular and new codon models of molecular evolution.

Documentation for [CodonMolecularEvolution](https://github.com/MurrellGroup/CodonMolecularEvolution.jl).

Descendant of [MolecularEvolution.jl](https://github.com/MurrellGroup/MolecularEvolution.jl), specializing in codon models.

### Collection of codon model methods
- [difFUBAR](@ref): Scalable Bayesian comparison of selection pressure
    - Perform a site-wise comparison of evolutionary pressure between two selected sets of branches.
    - Authors: Hassan Sadiq, Venkatesh Kumar, and Ben Murrell (original model development), Patrick Truong (benchmarking), Maximilian Danielsson (performance optimization).

### Design principles
- User-facing
    - Users with no Julia experience should be able to run these models.
- Scalability
    - Models should scale to large, real-world datasets. To keep the memory footprint down, we use `MolecularEvolution.jl`s `LazyPartition` when possible.
- Performance
    - We try to maintain competitive runtimes by using e.g. computational shortcuts and parallelization whenever we can.

### Package Authors and Maintainers
Maximilian Danielsson and Ben Murrell. Authors for specific models listed above.

```@index
```

```@autodocs
Modules = [CodonMolecularEvolution]
```

#=
CSV is included in the main module, but Plots is in an extension.
For each plotting function, foo, we:
    define foo(_, args...) as empty
    and in the extension, we'll define foo(::PlotBackend, args...)
=#

function plot_tagged_phylo_tree(args...; kwargs...) end

function FUBAR_plot_results(args...; kwargs...) end

function difFUBAR_plot_results(args...; kwargs...) end
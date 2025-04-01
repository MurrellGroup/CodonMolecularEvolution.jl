module FUBARShared
# This file should include all the stuff that is common across FUBAR methods.
# A FUBAR method at its most basic level should be able to take in a tagged tree
# and output an empirical prior distribution corresponding to a discretisation of the
# bivariate distribution over synonymous/non-synonymous substitution
struct FUBARgrid{T}
    grid_values::Vector{T}
    alpha_vec::Vector{T}
    beta_vec::Vector{T}
    alpha_ind_vec::Vector{Int64}
    beta_ind_vec::Vector{Int64}
    cond_lik_matrix::Matrix{T}
    LL_offset::T
    sites::Int64
    grid_function::Function
    LL_matrix::Matrix{T}
end

function gridplot(grid::FUBARgrid, θ; title = "")
    scatter(grid.alpha_ind_vec,grid.beta_ind_vec, zcolor = θ, c = :darktest, 
    markersize = sqrt(length(grid.alpha_ind_vec))/3.5, markershape=:square, markerstrokewidth=0.0, size=(400,350),
    label = :none, xticks = (1:length(grid.grid_values), round.(grid.grid_values,digits = 3)), xrotation = 90,
    yticks = (1:length(grid.grid_values), round.(grid.grid_values,digits = 3)), margin=6Plots.mm,
    xlabel = "α", ylabel = "β", title = title, colorbar = true, right_margin = 12Plots.mm)
    plot!(1:length(grid.grid_values),1:length(grid.grid_values),color = "grey", style = :dash, label = :none)
end

function FUBAR_init(treestring; verbosity=1, exports=true, disable_binarize=false, ladderize_tree = false)

    tree = gettreefromnewick(treestring, FelNode, disable_binarize=disable_binarize)
    if ladderize_tree
        MolecularEvolution.ladderize!(tree)
    end
    verbosity > 0 && println("Step 1: Initialization.")
    return tree
end


export 
    gridplot
end
using MolecularEvolution, CodonMolecularEvolution
using Test

@testset "CodonMolecularEvolution.jl" begin
    @testset "difFUBAR" begin
        include("difFUBAR_test.jl")
    end
end

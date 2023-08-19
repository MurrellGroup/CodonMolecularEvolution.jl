using MolecularEvolution, CodonMolecularEvolution
using Test

@testset "CodonMolecularEvolution.jl" begin
    # Write your tests here.
    @testset "difFUBAR" begin
        include("difFUBAR_test.jl")
    end
end

using Test
using FermionSimulation
using FermionSimulation: isevenperm

@testset "space config" begin
    sc = SpaceConfig((["↑", "↓"], ["1", "2", "3", "4"]))
    @test get_name_of_flavor(sc, 1) == ("↑", "1")
    @test get_name_of_flavor(sc, 2) == ("↓", "1")
    @test get_name_of_flavor(sc, 3) == ("↑", "2")
end

@testset "isevenperm" begin
    @test isevenperm([1, 2, 3, 4])
    @test !isevenperm([1, 3, 2, 4])
    @test isevenperm((1, 2, 3, 4))
    @test !isevenperm((1, 3, 2, 4))
end

@testset "fermionic operators" begin
    sc = SpaceConfig((["↑", "↓"], ["1", "2", "3", "4"]))
    @test FAnnilation(sc, 1) |> isfermionic
    @test FCreation(sc, 1) |> isfermionic
end

@testset "standard fermionic string" begin
    sc = SpaceConfig((["↑", "↓"], ["1", "2", "3", "4"]))
    @test_throws AssertionError StandardFermionicString(FAnnilation(sc, 1), FCreation(sc, 2)) isa StandardFermionicString{2, Tuple{FAnnilation{2}, FCreation{2}}}
    @test StandardFermionicString(FCreation(sc, 2), FAnnilation(sc, 1)) isa StandardFermionicString{2, Tuple{FCreation{2}, FAnnilation{2}}}
end

@testset "fermion hamiltonian" begin
    sc = SpaceConfig((["↑", "↓"], ["1", "2", "3", "4"]))
    fs1 = (FAnnilation(sc, 1), FCreation(sc, 2))
    fs2 = (FAnnilation(sc, 3), FCreation(sc, 4))
    h1 = FermionHamiltonian([1.0 => fs1])
    h2 = FermionHamiltonian([2.0 => fs2])
    @test h1 - h2 isa FermionHamiltonian
    @test push!(-h1, 3.0 => fs1) isa FermionHamiltonian
end

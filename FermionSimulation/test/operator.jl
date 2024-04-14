using Test
using FermionSimulation
using FermionSimulation: isevenperm

# @testset "space config" begin
#     sc = SpaceConfig((["↑", "↓"], ["1", "2", "3", "4"]))
#     @test get_name_of_flavor(sc, 1) == ("↑", "1")
#     @test get_name_of_flavor(sc, 2) == ("↓", "1")
#     @test get_name_of_flavor(sc, 3) == ("↑", "2")
# end

@testset "isevenperm" begin
    @test isevenperm([1, 2, 3, 4])
    @test !isevenperm([1, 3, 2, 4])
    @test isevenperm((1, 2, 3, 4))
    @test !isevenperm((1, 3, 2, 4))
end

@testset "fermionic operators" begin
    @test !(annilation(1).iscreation)
    @test creation(1).iscreation
end

@testset "standard fermionic string" begin
    fs = StandardFermionicString(1.0, annilation(1), creation(2))
    @test fs isa StandardFermionicString{Float64, 2}
    @test fs.coeff == -1.0
    @test fs.ops[1].flavor == 2
    @test fs.ops[2].flavor == 1
    fs2 = StandardFermionicString(2, creation(2), annilation(1))
    @test fs2 isa StandardFermionicString{Int, 2}
    @test fs2.coeff == 2
    @test fs2.ops[1].flavor == 2
    @test fs2.ops[2].flavor == 1
end

@testset "fermion hamiltonian" begin
    fs1 = StandardFermionicString(2.0, annilation(1), creation(2))
    fs2 = StandardFermionicString(3.0, annilation(3), creation(4))
    h1 = SimpleFermionHamiltonian([fs1])
    h2 = SimpleFermionHamiltonian([fs2])
    @test h1 - h2 isa SimpleFermionHamiltonian
    @test push!(-h1, fs1) isa SimpleFermionHamiltonian
end

@testset "upper_triangle_count" begin
    n=9
    for i in 1:(n-1)
        for j in (i+1):n
            i1,j1 = upper_triangle_count_inv(upper_triangle_count(i,j,n),n)
            @test i1 == i
            @test j1 == j
        end
    end
end

@testset "CSCFH quadratic" begin
    fs1 = StandardFermionicString(2.0, annilation(1), creation(2))
    fs2 = StandardFermionicString(3.0, annilation(3), creation(4))
    fs3 = StandardFermionicString(-1.0, creation(4), annilation(3))
    @show CSCFH([fs1,fs2,fs3],4)
end 

@testset "CSCFH quartic" begin
    fs1 = StandardFermionicString(2.0, annilation(4), creation(1), annilation(3), creation(4))
    fs2 = StandardFermionicString(3.0, annilation(3), creation(4), annilation(1), creation(2))
    fs3 = StandardFermionicString(-10.0, creation(4), annilation(3), creation(2), annilation(1))
    @show CSCFH([fs1,fs2,fs3],4)
end

@testset "CSCSimpleFermionHamiltonian" begin
    fs1 = StandardFermionicString(2.0, annilation(1), creation(2))
    fs2 = StandardFermionicString(3.0, annilation(3), creation(4))
    fs3 = StandardFermionicString(-1.0, creation(4), annilation(3))
    fs4 = StandardFermionicString(2.0, annilation(4), creation(1), annilation(3), creation(4))
    fs5 = StandardFermionicString(3.0, annilation(3), creation(4), annilation(1), creation(2))
    fs6 = StandardFermionicString(-10.0, creation(4), annilation(3), creation(2), annilation(1))
    @show CSCSimpleFermionHamiltonian(SimpleFermionHamiltonian([fs1,fs2,fs3],[fs4,fs5,fs6]),4)
end

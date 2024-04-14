using Test
using FermionSimulation

@testset "all_bases" begin
    bases = all_bases(3,5)
    @test bases[2] == ConstParticalNumberFermionBasis((1,2,4),)
    @test !(bases[2] == ConstParticalNumberFermionBasis((1,2,5),))
end

@testset "ConstParticalNumberFermionState" begin 
    psi = ConstParticalNumberFermionState([1.0 => ConstParticalNumberFermionBasis((1,2,3),), 2.0 => ConstParticalNumberFermionBasis((1,4,2),), 3.0 => ConstParticalNumberFermionBasis((1,2,5),)])
    @test -psi-psi isa ConstParticalNumberFermionState{Float64,3}
    @test 2.0 * psi isa ConstParticalNumberFermionState{Float64,3}
end

@testset "sort_and_combine_terms" begin
    psi = ConstParticalNumberFermionState([1.0 => ConstParticalNumberFermionBasis((1,6,3),), 2.0 => ConstParticalNumberFermionBasis((1,4,2),), 3.0 => ConstParticalNumberFermionBasis((1,2,5),), 4.0 => ConstParticalNumberFermionBasis((1,2,5),)])
    psi2 = ConstParticalNumberFermionState(sort_and_combine_terms(psi.megs))
    @test psi2.megs[2] == Pair(7.0, ConstParticalNumberFermionBasis((1,2,5),))
end

@testset "apply_SFS_on_basis" begin
    fs1 = StandardFermionicString(2.0, annilation(2), creation(1), annilation(7), creation(9), creation(2),annilation(9))
    base1 = ConstParticalNumberFermionBasis((2,3,7,9),)
    coeff,base2 = apply_SFS_on_basis(fs1,base1)
    @test coeff == 2.0
    @test base2 == ConstParticalNumberFermionBasis((1,2,3,9),)
end

@testset "apply_hamiltonian_quadratic" begin
    psi = ConstParticalNumberFermionState([1.0 => ConstParticalNumberFermionBasis((1,2,3),), 2.0 => ConstParticalNumberFermionBasis((1,2,4),), 3.0 => ConstParticalNumberFermionBasis((1,2,5),)])

    fs1 = StandardFermionicString(2.0, annilation(1), creation(2))
    fs2 = StandardFermionicString(3.0, annilation(3), creation(4))
    fs3 = StandardFermionicString(-1.0, creation(4), annilation(3))
    h=CSCFH([fs1,fs2,fs3],5)
    @show apply_hamiltonian(h, psi)
end


# h = (-2.0, (c6†, c1))+(-3.0, (c4†, c3))+(-1.0, (c3†, c4))+(2.0, (c3†, c6†, c1, c4))+(-3.0, (c3†, c6†, c1, c5))
#   = (-2.0, (c6†, c1))+(-4.0, (c4†, c3))+(2.0, (c3†, c6†, c1, c4))+(-3.0, (c3†, c6†, c1, c5))
# psi = 1.0|(1, 2, 3)> + 2.0|(1, 2, 4)> + 3.0|(1, 2, 5)>
# h * psi = -4.0|(1, 2, 4)> -7.0|(2, 3, 6)> -4.0|(2, 4, 6)> -6.0|(2, 5, 6)>
@testset "apply_hamiltonian" begin 
    psi = ConstParticalNumberFermionState([1.0 => ConstParticalNumberFermionBasis((1,2,3),), 2.0 => ConstParticalNumberFermionBasis((1,2,4),), 3.0 => ConstParticalNumberFermionBasis((1,2,5),)])

    fs1 = StandardFermionicString(2.0, annilation(1), creation(6))
    fs2 = StandardFermionicString(3.0, annilation(3), creation(4))
    fs3 = StandardFermionicString(-1.0, creation(4), annilation(3))
    fs4 = StandardFermionicString(2.0, annilation(4), creation(3), annilation(1), creation(6))
    fs5 = StandardFermionicString(3.0, annilation(1), creation(3), annilation(5), creation(6))
    h=CSCSimpleFermionHamiltonian(SimpleFermionHamiltonian([fs1,fs2,fs3],[fs4,fs5]),6)
    @test apply_hamiltonian(h, psi).megs[2] == Pair(-7.0, ConstParticalNumberFermionBasis((2,3,6),))
end

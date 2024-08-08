using Test
using FermionSimulation
using FermionSimulation.SparseArrays
using FermionSimulation: apply_SFS_on_basis_withoutsign, parity_judge

@testset "all_bases" begin
    bases = all_bases(3,5)
    @test bases[2] == FockBasis((1,2,4),)
    @test !(bases[2] == FockBasis((1,2,5),))
    dict = bases_dict(bases)
    @test dict[FockBasis((1,2,4),)] == 2
end

# h = (-2.0, (c5†, c1)) + (-3.0, (c4†, c3)) + (-1.0, (c4†, c3)) + (2.0, (c3†, c5†, c1, c4)) + (-3.0, (c3†, c5†, c1, c5))
#   = (-2.0, (c5†, c1)) + (-4.0, (c4†, c3)) + (2.0, (c3†, c5†, c1, c4)) + (-3.0, (c3†, c5†, c1, c5))
@testset "CSCSimpleFH_under_bases" begin
    fs1 = StandardFermionicString(2.0, annilation(1), creation(5))
    fs2 = StandardFermionicString(3.0, annilation(3), creation(4))
    fs3 = StandardFermionicString(-1.0, creation(4), annilation(3))
    fs4 = StandardFermionicString(2.0, annilation(4), creation(3), annilation(1), creation(5))
    fs5 = StandardFermionicString(-3.0, annilation(1), creation(3), creation(5), annilation(5))
    h=CSCSimpleFermionHamiltonian(SimpleFermionHamiltonian([fs1,fs2,fs3],[fs4,fs5]),5)
    bases = all_bases(3,5)
    dict = bases_dict(bases)
    @test CSCSimpleFH_under_bases(h, bases, dict) == sparse([8, 2, 9, 10, 6, 9, 8, 8, 10], [1, 1, 2, 4, 5, 8, 2, 3, 6], [-2.0, -4.0, -2.0, -2.0, -4.0, -4.0, 2.0, -3.0, 3.0], 10, 10)
end

@testset "apply_SFS_on_basis_withoutsign" begin
    fs = StandardFermionicString(2.0, annilation(2), creation(3), annilation(1), creation(5))
    basis = FockBasis((1,2,4),)
    flavors = ntuple(i->fs.ops[i].flavor, 4)
    @test apply_SFS_on_basis_withoutsign(flavors, basis) == FockBasis((5,3,4),)
end

@testset "parity_judge" begin
    fs = StandardFermionicString(2.0, annilation(2), creation(3), annilation(1), creation(5))
    fs3 = StandardFermionicString(-1.0, creation(4), annilation(3))
    basis = FockBasis((1,2,3),)
    flavors = ntuple(i->fs3.ops[i].flavor, 2)
    @test parity_judge(Float64,flavors, basis) == 1.0
end

@testset "kron_state" begin
    wave1 = [1.0, 2.0,0.0, 0.0, 0.0, 0.0]
    wave2 = [0.0, 0.0, 5.0, 6.0, 3.0, 4.0]
    wave3 = [0.0, 0.0, 0.0, 0.0, 3.0, 4.0]
    @test length(kron_state(wave1, wave2)) == 15
    @test length(kron_state(wave1, wave2, wave3)) == 20
end
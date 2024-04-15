using Test
using FermionSimulation
using FermionSimulation.Graphs
using FermionSimulation.LinearAlgebra

@testset "time_evolve" begin
    sg = SimpleGraph(Edge.([(1,2), (2,3), (3,4), (4,5), (5,6), (6,7)]))
    h = get_hamiltonian_from_sg(sg)
    n = nv(sg)
    hc = CSCSimpleFermionHamiltonian(h, n)
    bases = all_bases(1,n)
    dict = bases_dict(bases)
    H = CSCSimpleFH_under_bases(hc, bases, dict)
    random_wave =ones(ComplexF64, 7)
    random_wave = random_wave / norm(random_wave)
    waves = time_evolve(random_wave, ComplexF64.(H), 10, 1.0)
    @test waves[11] ==ComplexF64[-0.3320457865207849 + 0.01348357872068922im, -0.3599062873341397 + 0.026967157441378377im, -0.3877667881474948 - 0.09690113896565257im, -0.3877667881474948 - 0.220769435372684im, -0.3877667881474948 - 0.09690113896565256im, -0.3599062873341396 + 0.02696715744137864im, -0.33204578652078476 + 0.013483578720689197im]
end
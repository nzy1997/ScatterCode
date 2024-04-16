using Test
using FermionSimulation
using FermionSimulation.SparseArrays
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
    waves = []
    time_evolve(random_wave, ComplexF64.(H), 10, 1.0;cache=waves)
    @test waves[11] == ComplexF64[-0.3320457865207849 + 0.01348357872068922im, -0.3599062873341397 + 0.026967157441378377im, -0.3877667881474948 - 0.09690113896565257im, -0.3877667881474948 - 0.220769435372684im, -0.3877667881474948 - 0.09690113896565256im, -0.3599062873341396 + 0.02696715744137864im, -0.33204578652078476 + 0.013483578720689197im]
end

@testset "position_operators" begin
    bases = all_bases(2,5)
    dict = bases_dict(bases)
    @test position_operators(5, bases, dict) == SparseArrays.SparseMatrixCSC{Float64, Int64}[sparse([1, 2, 3, 4], [1, 2, 3, 4], [1.0, 1.0, 1.0, 1.0], 10, 10), sparse([1, 5, 6, 7], [1, 5, 6, 7], [1.0, 1.0, 1.0, 1.0], 10, 10), sparse([2, 5, 8, 9], [2, 5, 8, 9], [1.0, 1.0, 1.0, 1.0], 10, 10), sparse([3, 6, 8, 10], [3, 6, 8, 10], [1.0, 1.0, 1.0, 1.0], 10, 10), sparse([4, 7, 9, 10], [4, 7, 9, 10], [1.0, 1.0, 1.0, 1.0], 10, 10)]
end

@testset "expect_value" begin
    wave1 = [1.0, 2.0,0.0, 0.0, 0.0, 0.0]
    wave2 = [0.0, 0.0, 5.0, 6.0, 0.0, 0.0]
    wave3 = [0.0, 0.0, 0.0, 0.0, 3.0, 4.0]
    psi = kron_state(wave1, wave2, wave3)

    bases = all_bases(3,6)
    dict = bases_dict(bases)
    operators = position_operators(6, bases, dict)
    @show [expect_value(psi, operator) for operator in operators]
    @test expect_value(psi, operators[1]) == (15^2+20^2+18^2+24^2)
end
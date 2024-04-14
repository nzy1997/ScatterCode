using Test
using FermionSimulation
using FermionSimulation.Graphs

@testset "get_hamiltonian_from_sg" begin
    sg = SimpleGraph(Edge.([(1,3), (2,3), (3,4), (4,5), (5,6), (6,7), (7,4)]))
    h = get_hamiltonian_from_sg(sg)
    hc=CSCSimpleFermionHamiltonian(h, 7)
    bases = all_bases(1,7)
    dict = bases_dict(bases)
    @test CSCSimpleFH_under_bases(hc, bases, dict) == -adjacency_matrix(sg)/2
end
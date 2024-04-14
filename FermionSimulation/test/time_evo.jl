using Test
using FermionSimulatio
using FermionSimulation.Graphs
using FermionSimulation.LinearAlgebra

@testset "time_evo" begin
    sg = SimpleGraph(Edge.([(1,2), (2,3), (3,4), (4,5), (5,6), (6,7)]))
    h = get_hamiltonian_from_sg(sg)
    n = nv(sg)
    hc = CSCSimpleFermionHamiltonian(h, n)
    bases = all_bases(1,n)
    dict = bases_dict(bases)
    H = CSCSimpleFH_under_bases(hc, bases, dict)
    random_wave = rand(ComplexF64, n)
    random_wave = random_wave / norm(random_wave)
    waves = time_evo(random_wave, ComplexF64.(H), 10, 1.0)
    @test waves[11] == ComplexF64[-0.34905597635463975 + 0.06710690636041369im, -0.26288076187405673 - 0.20358592003901935im, -0.040195488374988306 - 0.21029204557207012im, -0.05096789256656367 - 0.4664248977957689im, -0.34762555070892703 - 0.44379597926331443im, -0.16084710894666532 - 0.3309250104876749im, -0.18471295232409812 - 0.09906719497463454im]
end
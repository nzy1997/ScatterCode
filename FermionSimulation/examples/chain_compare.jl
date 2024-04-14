using Scattering
using Scattering.Graphs
using Scattering.LinearAlgebra

gt = graph_with_tails(SimpleGraph(Edge.([(1,2)])); heads=[2], tailsize=499)
h0 = Diagonal(zeros(2))
waves = Vector{ComplexF64}[]
v = simulate_graph_with_tails(gt, h0; Nt=1000, cache=waves)
animate_wave(gt, waves, step = 100, pathname = "examples/chain1.gif")

using FermionSimulation
h = get_hamiltonian_from_sg(gt.graph)
n = nv(gt.graph)
hc=CSCSimpleFermionHamiltonian(h, n)
bases = all_bases(1,n)
dict = bases_dict(bases)
H = CSCSimpleFH_under_bases(hc, bases, dict)
waves2 = time_evo(waves[1], ComplexF64.(H), 998, 1.0)
animate_wave(gt, waves2, step = 100, pathname = "examples/chain2.gif")
isapprox(waves,waves2;atol=1e-8)
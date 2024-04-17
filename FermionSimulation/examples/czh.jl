using FermionSimulation
using FermionSimulation.Scattering
using FermionSimulation.Graphs

k0 = - π/2
gt,n1,n2 = example_gt_cz(200,50)
hs=SimpleFermionHamiltonian([StandardFermionicString(-5.0, creation(n1+i),creation(n2+i),annilation(n1+i), annilation(n2+i)) for i in -25:25])
waves = simulate_gt(gt, k0, [1,3];Nt = 499, tag = true,hs)

[add_edge!(gt.graph, n1+i, n2+i) for i in -25:25]
animate_graph(gt,waves, step = 25, pathname = "examples/chaincz3.gif")

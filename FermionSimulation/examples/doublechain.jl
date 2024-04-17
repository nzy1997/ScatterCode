using FermionSimulation
using FermionSimulation.Scattering
using FermionSimulation.Graphs

k0 = - π/3
n1 = 20

gt = example_gt_doublechain(n1)
waves0 = simulate_gt(gt, k0, [1,2];Nt = 490, tag = false)
hs=SimpleFermionHamiltonian([StandardFermionicString(-1000, creation(n1+1),creation(1),annilation(n1+1), annilation(1))])
waves = simulate_gt(gt, k0, [1,2];Nt = 490, tag = false,hs)

add_edge!(gt.graph, 1, n1+1)
animate_graph(gt,waves, step = 25, pathname = "examples/chaindb2.gif")

using FermionSimulation
using FermionSimulation.Scattering
using FermionSimulation.Graphs

gt = graph_with_tails(SimpleGraph(Edge.([(1,2)])); heads=[2], tailsize=100)

hs=SimpleFermionHamiltonian([StandardFermionicString(-100.0, creation(250), annilation(250))])
waves = simulate_gt(gt, k0, [1];Nt = 100, tag = false)

animate_wave(gt,waves, step = 5, pathname = "examples/chain.gif")

using FermionSimulation
using FermionSimulation.Scattering
using FermionSimulation.Graphs

gt = graph_with_tails(SimpleGraph(Edge.([(1,2)])); heads=[2], tailsize=498)

hs=SimpleFermionHamiltonian([StandardFermionicString(-100.0, creation(250), annilation(250))])
waves = simulate_gt(gt, k0, [1];Nt = 499, tag = false,hs)

animate_wave(gt,waves, step = 25, pathname = "examples/reflect.gif")

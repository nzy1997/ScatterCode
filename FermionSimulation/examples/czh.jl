using FermionSimulation
using FermionSimulation.Scattering
using FermionSimulation.Graphs

# k0 = - π/2
# gt,n1,n2 = example_gt_cz(200,50)
# hs=SimpleFermionHamiltonian([StandardFermionicString(-5.0, creation(n1+i),creation(n2+i),annilation(n1+i), annilation(n2+i)) for i in -25:25])
# [add_edge!(gt.graph, n1+i, n2+i) for i in -25:25]



gt = example_gt_mss(100,100)
waves = simulate_gt(gt, [-pi/2,-pi/4], [1,3]; Nt = 499, tag = true)

animate_wave(gt,waves, step = 25, pathname = "examples/chainmss.gif")

animate_graph(gt,waves, step = 25, pathname = "examples/chainmss4.gif")

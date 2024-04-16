using FermionSimulation
using FermionSimulation.Scattering
using FermionSimulation.Graphs

k0 = - π/2
gt,n1,n2 = example_gt_cz(200,50)
fs=StandardFermionicString(-1.0, creation(n1),creation(n2),annilation(n1), annilation(n2))
waves = simulate_gt(gt, k0, [1,3];Nt = 499, tag = true,fs)

animate_graph(gt, waves, step = 25, pathname = "examples/chaincz.gif")

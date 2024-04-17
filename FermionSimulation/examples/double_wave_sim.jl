using FermionSimulation
using FermionSimulation.Scattering
using FermionSimulation.Graphs
using FermionSimulation.Scattering: gaussian_packet_ontail

gt = graph_with_tails(SimpleGraph(Edge.([(1,2)])); heads = [1,2], tailsize = 20)
k0 = π/4
waves = simulate_gt(gt, k0, [1,2];tag = true)

animate_wave(gt, 5*waves, step = 100, pathname = "examples/chain2pn.gif", framerate=5)
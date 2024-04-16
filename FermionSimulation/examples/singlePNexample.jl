using FermionSimulation

using FermionSimulation.Scattering
using FermionSimulation.Graphs
using FermionSimulation.Scattering: gaussian_packet_ontail

gt = example_gt_h(200,50)
k0 = - π/2
waves = simulate_gt(gt, k0, [1];Nt = 499, tag = false)

animate_wave(gt, waves, step = 25, pathname = "examples/chain2.gif")


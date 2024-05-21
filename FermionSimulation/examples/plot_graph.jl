using FermionSimulation.Graphs
using Plots, GraphRecipes
using FermionSimulation
using FermionSimulation: momentum_switch_sg!,add_chain!

g = example_gt_cz(2,4).graph
g = example_gt_doublechain(4).graph

g= momentum_switch_sg!()
add_chain!(g,1,2,10)

g =example_gt_mss(4,6).graph


graphplot(g, curves=false)

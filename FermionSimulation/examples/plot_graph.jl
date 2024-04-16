using FermionSimulation.Graphs
using Plots, GraphRecipes
using FermionSimulation

g = example_gt_cz(2,4).graph
graphplot(g, curves=false)

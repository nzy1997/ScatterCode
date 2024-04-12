using Test, Scattering, Scattering.Graphs, LinearAlgebra
using Scattering.LuxorGraphPlot: Luxor


g = andrew_momentum_separator()
waves, gt = simulate_ScatterGraph(g, -3/4;tailsize=499)
animate_wave(gt, waves, step = 100, pathname = "examples/chain-0.75.gif")

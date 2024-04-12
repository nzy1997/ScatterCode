using Scattering
using Scattering.Graphs
using LinearAlgebra
using Scattering.LuxorGraphPlot: Luxor
using Random

Random.seed!(435)
x0 = rand(Float64, 20)
g = andrew_hadamard(x0)
outcome = optimize_weighted_momentum(g,rand(Float64))

waves,gt=simulate_ScatterGraph(outcome.g, outcome.minimizer[1])
animate_wave(gt, waves, step=100,pathname="examples/chain.gif") isa Luxor.AnimatedGif

gh = andrew_hadamard(ones(20))
waves,gt=simulate_ScatterGraph(gh, -0.5)
animate_wave(gt, waves, step=100,pathname="examples/chain.gif") isa Luxor.AnimatedGif

optimize_momentum(outcome.g,-0.7)
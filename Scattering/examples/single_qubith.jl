include("check_graph.jl")

Random.seed!(989117422861198)
x0 = rand(Float64, 20)
g = andrew_hadamard(x0)
outcome = optimize_weighted_momentum(g, rand(Float64))

waves, gt = simulate_ScatterGraph(outcome.g, -outcome.minimizer[1])
animate_wave(gt, waves, step = 100, pathname = "examples/chain.gif") isa Luxor.AnimatedGif

ulist= [get_u1(outcome.s), get_u2(outcome.s)]
universal_check(ulist,1)

gms = andrew_momentum_separator()
plot_rr_momentum(outcome.g, "examples/rr_momentum.png")


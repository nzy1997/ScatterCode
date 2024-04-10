using Scattering
using Scattering.Graphs
using LinearAlgebra
using Scattering.LuxorGraphPlot: Luxor

# sg=SimpleGraph(Edge.([(1,5), (2,7), (3, 8), (4, 10), (5,6), (5,8), (6,7), (7,10), (8,9), (9, 10)]))
sg = SimpleGraph(Edge.([(1,7), (5,6),(4,6),(3,11),(6,7),(6,10),(5,9),(5,8),(7,8),(7,10),(8,9),(9,10),(8,11),(9,13),(10,11),(10,13),(8,12),(11,12),(2,12),(12,13)]))

g = ScatterGraph(sg, [1, 2], [3, 4])
x0 = rand(Float64, length(edges(g.graph))+1)
x, min = optimize_weighted_momentum(g, x0)

z=exp(pi*im*x[1])
s=scatter_matrix(g, z, x[2:end])
isapprox(s[1:2,1:2], zeros(2,2), atol=1e-5)
isapprox(s[3:4,3:4], zeros(2,2), atol=1e-5)
u1=s[1:2,3:4]
u2=s[3:4,1:2]
ulist = [u1, u2]


gt = graph_with_tails(sg; heads=[1, 2,3,4], tailsize=499)
W = matrix_adjacency(g.graph, x[2:end])
h0 = W-matrix_adjacency(g.graph, ones(Float64, ne(g.graph)))
h0=h0/2
# h0 = Diagonal(zeros(13))
waves = Vector{ComplexF64}[]
v = simulate_graph_with_tails(gt, h0;k0=x[1]*pi, Nt=1500, cache=waves,intail=2)
animate_wave(gt, waves, step=100) isa Luxor.AnimatedGif
plot_transmission(andrew_momentum_separator(); outtail=4) isa Luxor.Drawing
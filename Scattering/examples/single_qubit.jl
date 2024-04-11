using Scattering
using Scattering.Graphs
using LinearAlgebra
using Scattering.LuxorGraphPlot: Luxor
using Random

struct OptimizedResult{RT}
    minimizer::Vector{RT}
    minimum::RT
    s::Matrix{Complex{RT}}
    g::ScatterGraph
end
function get_u1(r::OptimizedResult)
    n = size(r.s, 1) ÷ 2
    return r.s[1:n,n+1:end]
end
function get_u2(r::OptimizedResult)
    n = size(r.s, 1) ÷ 2
    return r.s[n+1:end,1:n]
end
function get_r1(r::OptimizedResult)
    n = size(r.s, 1) ÷ 2
    return r.s[1:n,1:n]
end
relection_rate(r::OptimizedResult) = sum(abs2, get_r1(r))

function generate_model_optim()
    Random.seed!(1234)
    sg = SimpleGraph(Edge.([(1,7), (5,6),(4,6),(3,11),(6,7),(6,10),(5,9),(5,8),(7,8),(7,10),(8,9),(9,10),(8,11),(9,13),(10,11),(10,13),(8,12),(11,12),(2,12),(12,13)]))
    x0 = rand(Float64, length(edges(g.graph)))
    inputs, outputs = [1, 2], [3, 4]
    return ScatterGraph(sg, x0, inputs, outputs)
end

function optimize_weighted_momentum(sg::SimpleGraph, weights::AbstractVector, inputs, outputs)
    g = ScatterGraph(sg, inputs, outputs)
    x, min = optimize_weighted_momentum(g, weights)
    z = exp(pi*im*x[1])
    s = scatter_matrix(g, z, x[2:end])
    return OptimizedResult(x, min, s, g)
end

x,min,u1,u2,s,sg,g=generate_model_optim()
gt = graph_with_tails(sg; heads=[1, 2,3,4], tailsize=500)

function simu(gt,x,g)
    h0 = matrix_adjacency(g.graph, x[2:end] .- 1) ./ 2  # difference in weight

    H = get_hamiltonian(gt, h0)
    any(==(1), H[gt.center, gt.center]) && error("Center has self-loop")

    k0 = (mod(x[1],2.0) > 1.0 ? (mod(x[1],2.0)-2.0)*pi : mod(x[1],2.0)*pi)

    # h0 = Diagonal(zeros(13))
    # k0 =-pi/2

    waves = Vector{ComplexF64}[]
    simulate_graph_with_tails(gt, h0;k0, Nt=1500, cache=waves,intail=1)
    return waves
end
# x=[-1/2,ones(20)...]

waves=simu(gt,x,g)
animate_wave(gt, waves, step=100,pathname="examples/chain.gif") isa Luxor.AnimatedGif
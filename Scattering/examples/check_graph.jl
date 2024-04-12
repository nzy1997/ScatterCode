using Scattering
using Scattering.Graphs
using LinearAlgebra
using Scattering.LuxorGraphPlot: Luxor
using Random
using CairoMakie

function plot_rr_momentum(g::ScatterGraph, save_name::String)
	k = -1+1e-10:0.003:-1e-10
	# k=1e-10:0.003:1-1e-10
	rr = [relection_rate(scatter_matrix(g, exp(im * ki * pi))) for ki in k]
	fig = Figure()
	ax = Axis(fig[1, 1])
	lines!(ax, k, rr)
	save(save_name, fig)
	return fig
end

function check_sg(sg::SimpleGraph,input_vertex::Vector{Int}, output_vertex::Vector{Int},save_path::String, iters::Int=10)
    for i in 1:iters
        s = rand(1:10^15)
        Random.seed!(s)
        x0 = rand(Float64, ne(sg))
        g = ScatterGraph(sg, x0, input_vertex, output_vertex)
        x1=rand(Float64)
        outcome = optimize_weighted_momentum(g, x1)
        @show s outcome.minimum
        if outcome.minimum < 1e-10
            waves, gt = simulate_ScatterGraph(outcome.g, outcome.minimizer[1])
            animate_wave(gt, waves, step = 100, pathname = joinpath(save_path,"$s.gif")) isa Luxor.AnimatedGif
            plot_rr_momentum(outcome.g, joinpath(save_path,"$s.png"))
        end
    end
end

# g = andrew_hadamard()
# s = 756080199614591, 3.839080041220371e-15
# s = 356850098191873, 6.635008359660213e-14
# check_sg(g.graph, g.input_vertex, g.output_vertex, "examples/andrew_hadamard/x1",100)
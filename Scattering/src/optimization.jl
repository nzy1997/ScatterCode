struct OptimizedResult{RT}
    minimizer::Vector{RT}
    minimum::RT
    s::Matrix{Complex{RT}}
    g::ScatterGraph
end

function train_loss(g::ScatterGraph)
	return x -> abs_loss(g, exp(pi * im * x[1]))
end

function train_weighted_loss(g::ScatterGraph)
	return x -> abs_loss(g.graph, x[2:end], g.input_vertex, g.output_vertex, exp(pi * im * x[1]))
end

function weighted_gra!(G, x0, g::ScatterGraph)
	G[1:length(x0)] = ForwardDiff.gradient(train_weighted_loss(g), x0)
	return G
end

function grad(g::ScatterGraph)
	return (G, x) -> weighted_gra!(G, x, g)
end

function optimize_momentum(g::ScatterGraph, x0::Real)
	a = optimize(train_loss(g), [x0])
	return a.minimizer, a.minimum
end

function optimize_weighted_momentum(g::ScatterGraph, k0)
	a = optimize(train_weighted_loss(g), grad(g), [k0,g.weights...], LBFGS())
	g2 = ScatterGraph(g.graph, a.minimizer[2:end], g.input_vertex, g.output_vertex) 
	return OptimizedResult(a.minimizer, a.minimum, scatter_matrix(g2, exp(pi * im * a.minimizer[1])), g2)
end
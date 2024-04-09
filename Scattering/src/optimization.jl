function train_loss(graph::ScatterGraph)
	return x -> abs2_loss(graph, exp(pi * im * x[1]))
end

function train_weighted_loss(g::ScatterGraph)
	return x -> abs2_loss(g, exp(pi * im * x[1]), x[2:end])
end

function weighted_gra!(G, x0, g::ScatterGraph)
	G[1:length(x0)] = ForwardDiff.gradient(train_weighted_loss(g), x0)
end

function grad(g::ScatterGraph)
	return (G, x) -> weighted_gra!(G, x, g)
end

function optimize_momentum(g::ScatterGraph)
	# a = optimize(train_loss(g), grad(g), [0.0], LBFGS())
	a = optimize(train_loss(g), [0.0])
	return a.minimizer, a.minimum
end

function optimize_weighted_momentum(g::ScatterGraph, x0)
	a = optimize(train_weighted_loss(g), grad(g), x0, LBFGS())
	# a = optimize(train_weighted_loss(g), ones(Float64, length(edges(g.graph))+1))  
	return a.minimizer, a.minimum
end
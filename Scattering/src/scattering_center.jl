struct ScatterGraph{RT}
	graph::SimpleGraph{Int}
	weights::Vector{RT}
	input_vertex::Vector{Int}
	output_vertex::Vector{Int}
	function ScatterGraph(graph::SimpleGraph{Int}, weights::Vector{RT}, input_vertex::Vector{Int}, output_vertex::Vector{Int}) where RT
		@assert length(weights) == ne(graph)
		@assert all(1 <= i <= nv(graph) for i in input_vertex)
		@assert all(1 <= i <= nv(graph) for i in output_vertex)
		@assert input_vertex ∩ output_vertex |> isempty
		new(graph, weights, input_vertex, output_vertex)
	end
end

function matrix_adjacency(g::SimpleGraph{Int}, weights::Vector{T}) where T
	@assert length(weights) == ne(g)
	edge = collect(edges(g))
    W = zeros(T,(nv(g),nv(g)))
	@inbounds for i in 1:length(weights)
		W[edge[i].src, edge[i].dst] = weights[i]
		W[edge[i].dst, edge[i].src] = weights[i]
	end
	return W
end

function create_Q(nv::Int, input_vertices, output_vertices)
	Q = Diagonal(ones(Int, nv))
	for i in input_vertices
		Q.diag[i] = 0
	end
	for i in output_vertices
		Q.diag[i] = 0
	end
	return Q
end

function abs_loss(g::ScatterGraph, z::T, W::Matrix) where T
	nv = size(W, 1)
	Q = create_Q(nv, g.input_vertex, g.output_vertex)
	A = I - W * z + z^2 * Q
	q, r = qr(A)
	loss = 0.0
	for k in 1:length(g.input_vertex)
		b = zeros(T, nv)
		b[g.input_vertex[k]] = 1 - z^2
		x = UpperTriangular(r) \ (q' * b)
		loss += sum(abs2.(x[g.input_vertex]))
		loss -= abs2(x[g.input_vertex[k]])
		loss += abs2(x[g.input_vertex[k]] - 1)
	end
	return loss
end

function abs_loss(g::ScatterGraph, z)
	W = matrix_adjacency(g.graph, g.weights)
    return abs_loss(g, z, W)
end

# Scatter Matrix
function perm_adjacency(g::ScatterGraph, W::Matrix)
	perm = [g.input_vertex..., g.output_vertex..., setdiff(1:size(W, 1), vcat(g.input_vertex, g.output_vertex))...]
	return W[perm, perm]
end

Q(A, B, D, z) = I - z * (A + B' * inv(I * (1 / z + z) - D) * B)
function scatter_matrix(g::ScatterGraph, z, W::Matrix)
	N = length(g.input_vertex) + length(g.output_vertex)
	A = W[1:N, 1:N]
	B = W[(N+1):end, 1:N]
	D = W[(N+1):end, (N+1):end]
	return -inv(Q(A, B, D, z)) * Q(A, B, D, inv(z))
end

function scatter_matrix(g::ScatterGraph, z)
	return scatter_matrix(g, z, ones(Float64, ne(g.graph)))
end

function scatter_matrix(g::ScatterGraph, z, weights::Vector)
	W = matrix_adjacency(g.graph, weights)
    return scatter_matrix(g, z, perm_adjacency(g, W))
end
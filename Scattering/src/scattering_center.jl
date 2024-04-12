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
		new{RT}(graph, weights, input_vertex, output_vertex)
	end
end
ScatterGraph(graph::SimpleGraph{Int},  input_vertex::Vector{Int}, output_vertex::Vector{Int})=ScatterGraph(graph, ones(Float64, ne(graph)), input_vertex, output_vertex)

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
function abs_loss(g::ScatterGraph,z)
	return abs_loss(g.graph, g.weights, g.input_vertex, g.output_vertex, z)
end

function abs_loss(graph::SimpleGraph, weights, input_vertex, output_vertex, z::T) where T
	W = matrix_adjacency(graph, weights)
	nv = size(W, 1)
	Q = create_Q(nv, input_vertex, output_vertex)
	A = I - W * z + z^2 * Q
	q, r = qr(A)
	loss = 0.0
	for k in input_vertex
		b = zeros(T, nv)
		b[k] = 1 - z^2
		x = UpperTriangular(r) \ (q' * b)
		loss += sum(abs.(x[input_vertex]))
		loss -= abs(x[k])
		loss += abs(x[k] - 1)
	end
	return loss
end

# Scatter Matrix
function perm_adjacency(g::ScatterGraph, W::Matrix)
	perm = [g.input_vertex..., g.output_vertex..., setdiff(1:size(W, 1), vcat(g.input_vertex, g.output_vertex))...]
	return W[perm, perm]
end

Q(A, B, D, z) = I - z * (A + B' * inv(I * (1 / z + z) - D) * B)
function scatter_matrix(g::ScatterGraph, z)
	W0 = matrix_adjacency(g.graph, g.weights)
	perm = [g.input_vertex..., g.output_vertex..., setdiff(1:size(W0, 1), vcat(g.input_vertex, g.output_vertex))...]
	W = W0[perm, perm]
	N = length(g.input_vertex) + length(g.output_vertex)
	A = W[1:N, 1:N]
	B = W[(N+1):end, 1:N]
	D = W[(N+1):end, (N+1):end]
	Q1 = Q(A, B, D, z)
	Q2 = Q(A, B, D, inv(z))
	# c1, c2 = cond(Q1), cond(Q2)
	# (c1 > 10 || c2 > 10) && throw(ArgumentError("Condition number of Q1 and Q2 should be less than 10, got $c1 and $c2"))
	# @show c1, c2,cond(I * (1 / z + z) - D)
	return -inv(Q1) * Q2
end

function get_u1(s::AbstractMatrix)
    n = size(s, 1) ÷ 2
    return s[1:n,n+1:end]
end

function get_u2(s::AbstractMatrix)
    n = size(s, 1) ÷ 2
    return s[n+1:end,1:n]
end

function get_r1(s::AbstractMatrix)
    n = size(s, 1) ÷ 2
    return s[1:n,1:n]
end

relection_rate(s::AbstractMatrix) = sum(abs2, get_r1(s))
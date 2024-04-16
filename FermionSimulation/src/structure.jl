function get_hamiltonian_from_sg(sg::SimpleGraph)
	quadratic = StandardFermionicString{Float64, 2}[]
	for edge in edges(sg)
		push!(quadratic, StandardFermionicString(0.5, creation(edge.dst), annilation(edge.src)))
		push!(quadratic, StandardFermionicString(0.5, creation(edge.src), annilation(edge.dst)))
	end
	return SimpleFermionHamiltonian(quadratic)
end

function get_CSCSimpleFH(sg::SimpleGraph, PN::Int)
	h = get_hamiltonian_from_sg(sg)
	n = nv(sg)
	hc = CSCSimpleFermionHamiltonian(h, n)
	bases = all_bases(PN, n)
	dict = bases_dict(bases)
	return CSCSimpleFH_under_bases(hc, bases, dict), bases, dict, n
end

function andrew_hadamard_sg!(startpoint::Int = 1, sg::SimpleGraph = SimpleGraph(13))
	n = startpoint - 1
	add_edge!(sg, n + 1, n + 7)
	add_edge!(sg, n + 5, n + 6)
	add_edge!(sg, n + 4, n + 6)
	add_edge!(sg, n + 3, n + 11)
	add_edge!(sg, n + 6, n + 7)
	add_edge!(sg, n + 6, n + 10)
	add_edge!(sg, n + 5, n + 9)
	add_edge!(sg, n + 5, n + 8)
	add_edge!(sg, n + 7, n + 8)
	add_edge!(sg, n + 7, n + 10)
	add_edge!(sg, n + 8, n + 9)
	add_edge!(sg, n + 9, n + 10)
	add_edge!(sg, n + 8, n + 11)
	add_edge!(sg, n + 9, n + 13)
	add_edge!(sg, n + 10, n + 11)
	add_edge!(sg, n + 10, n + 13)
	add_edge!(sg, n + 8, n + 12)
	add_edge!(sg, n + 11, n + 12)
	add_edge!(sg, n + 2, n + 12)
	add_edge!(sg, n + 12, n + 13)
	return sg
end

function example_graph_h(tailsize::Int, middle_tail::Int)
	sg = andrew_hadamard_sg!()
	add_vertices!(sg, 13)
	sg = andrew_hadamard_sg!(14, sg)
	sg2 = graph_with_tails(sg ; heads = [3,4], tailsize= 1+2* middle_tail).graph
    sg2 = graph_with_tails(sg2 ; heads = [16,17], tailsize).graph
	add_edge!(sg2, 26 + 1+2*middle_tail, 14)
	add_edge!(sg2, 26 + 2+4 * middle_tail, 15)
	return sg2
end
example_gt_h(tailsize::Int, middle_tail::Int) = graph_with_tails(example_graph_h(tailsize, middle_tail); heads = [1,2], tailsize)

function example_gt_cz(tailsize::Int, middle_tail::Int)
    sg = example_graph_h(tailsize, middle_tail)
    add_vertices!(sg, 1)
    n = nv(sg)
    sg2 = graph_with_tails(sg ; heads = [n], tailsize=1+2*middle_tail+tailsize + 5).graph
    # add_edge!(sg2, 26 + 3 * middle_tail+2, n+2+middle_tail+1)
    return graph_with_tails(sg2; heads = [1,2,n], tailsize), 26 + 3 * middle_tail+2, n+2+middle_tail+1
end
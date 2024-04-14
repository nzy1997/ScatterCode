function get_hamiltonian_from_sg(sg::SimpleGraph)
    quadratic = StandardFermionicString{Float64, 2}[]
    for edge in edges(sg)
        push!(quadratic, StandardFermionicString(0.5, creation(edge.dst), annilation(edge.src)))
        push!(quadratic, StandardFermionicString(0.5, creation(edge.src), annilation(edge.dst)))
    end
    return SimpleFermionHamiltonian(quadratic)
end
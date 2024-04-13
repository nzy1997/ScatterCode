struct FermionHamiltonian{T}
    qubit_number::Int
    two_body_terms::Vector{Tuple{Int,Int,T}}
    four_body_terms::Vector{Tuple{Int,Int,Int,Int,T}}
end

function get_hamiltonian(fh::FermionHamiltonian,basis::Vector{FermionBasis{n,N}}) where {n,N}
    H = zeros(Float64,length(basis),length(basis))

    return 
end

function get_hamiltonian(fh::FermionHamiltonian,n::Int)
    bases = all_bases(n,fh.qubit_number)
    return get_hamiltonian(fh, bases)
end

function get_hamiltonian(sg::SimpleGraph,n::Int)
    bases = all_bases(n,nv(sg))
    return get_hamiltonian(sg, bases)
end

function get_hamiltonian(sg::SimpleGraph,basis::Vector{FermionBasis{n,N}}) where {n,N}
    
    return n,N
end
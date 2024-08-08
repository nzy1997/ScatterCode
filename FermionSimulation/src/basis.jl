abstract type AbstractBasis end
abstract type AbstractFermionBasis <: AbstractBasis end
abstract type AbstractFermionState end

struct ConstParticalNumberFermionBasis{PN} <: AbstractFermionBasis
    indices::NTuple{PN,Int}
    function ConstParticalNumberFermionBasis(indices::NTuple{PN,Int}) where PN
        @assert length(indices) == PN "The number of indices should be same as the number of particals"
        new{PN}(TupleTools.sort(indices))
    end
end
Base.:(==)(b1::ConstParticalNumberFermionBasis{PN}, b2::ConstParticalNumberFermionBasis{PN}) where PN = all(b1.indices == b2.indices)

function Base.show(io::IO,::MIME"text/plain", c::ConstParticalNumberFermionBasis{PN}) where PN
    print(io, "|$(c.indices)>")
end
Base.show(io::IO, c::ConstParticalNumberFermionBasis) = Base.show(io, MIME("text/plain"), c)

function all_bases(PN::Int,n::Int)
    bases = ConstParticalNumberFermionBasis{PN}[]
    for pos in combinations(1:n, PN)
        push!(bases, ConstParticalNumberFermionBasis((pos...,)))
    end
    return bases
end

struct ConstParticalNumberFermionState{T,PN} <: AbstractFermionState
    megs::Vector{Pair{T,ConstParticalNumberFermionBasis{PN}}}
    function ConstParticalNumberFermionState(megs::Vector{Pair{T,ConstParticalNumberFermionBasis{PN}}}) where {T,PN}
        new{T,PN}(sort_and_combine_terms(megs))
    end
end

function Base.show(io::IO,::MIME"text/plain", s::ConstParticalNumberFermionState{T,PN}) where {T,PN}
    for (coeff, basis) in s.megs
        println(io, coeff, basis)
    end
end
Base.show(io::IO,s::ConstParticalNumberFermionState) = Base.show(io, MIME("text/plain"), s)

Base.:+(s1::ConstParticalNumberFermionState{T,PN}, s2::ConstParticalNumberFermionState{T,PN}) where {T,PN} = ConstParticalNumberFermionState(vcat(s1.megs, s2.megs))
Base.:*(x::Number, s::ConstParticalNumberFermionState{T,PN}) where {T,PN} = ConstParticalNumberFermionState([Pair(x * coeff, basis) for (coeff, basis) in s.megs])
Base.:-(s::ConstParticalNumberFermionState{T,PN}) where {T,PN} = (-1) * s
Base.:-(s1::ConstParticalNumberFermionState{T,PN}, s2::ConstParticalNumberFermionState{T,PN}) where {T,PN} = s1 + (-s2)

function basis_sort(bs1::ConstParticalNumberFermionBasis{PN},bs2::ConstParticalNumberFermionBasis{PN}) where PN
    for i in 1:PN
        bs1.indices[i] < bs2.indices[i] && return true
        bs1.indices[i] > bs2.indices[i] && return false
    end
    return false
end

function sort_and_combine_terms(megs::Vector{Pair{T,ConstParticalNumberFermionBasis{PN}}}) where {T,PN}
    order = sortperm(1:length(megs), lt = (i,j) -> basis_sort(megs[i][2],megs[j][2]))
    new_megs = Pair{T,ConstParticalNumberFermionBasis{PN}}[]
    base_pre = ConstParticalNumberFermionBasis((fill(1,PN)...,))
    for index in order
        coeff, basis = megs[index]
        if basis == base_pre
            c = new_megs[end].first + coeff
            pop!(new_megs)
            push!(new_megs, Pair(c, basis))
        else
            push!(new_megs, Pair(coeff, basis))
            base_pre = basis
        end
    end
    return new_megs
end
sort_and_combine_terms(psi::ConstParticalNumberFermionState) = ConstParticalNumberFermionState(sort_and_combine_terms(psi.megs))

function apply_SFS_on_basis(fs::StandardFermionicString{T,N}, basis::ConstParticalNumberFermionBasis) where {T,N}
    coeff = fs.coeff
    for op in fs.ops[N:-1:(N÷2+1)]
        index=findfirst(x -> x == op.flavor, basis.indices) 
        isnothing(index) && return zero(T), basis
        coeff *= (-1)^(index+1)
    end
    new_basis = basis.indices[findall(x -> x ∉ [fs.ops[i].flavor for i in N:-1:(N÷2+1)], basis.indices)]
    for op in fs.ops[1:N÷2]
        index = findfirst(x -> x >= op.flavor, new_basis) 
        if isnothing(index)
            coeff *= (-1)^length(new_basis)
        else
            isequal(new_basis[index],op.flavor) && return zero(T), basis    
            coeff *= (-1)^(index+1)
        end
    end
    new_basis=(sort(union([fs.ops[i].flavor for i in 1:N÷2], new_basis))...,)
    return coeff, ConstParticalNumberFermionBasis(new_basis)
end

function apply_hamiltonian(h::CSCFH{T,N}, state::ConstParticalNumberFermionState{T,PN}) where {N,T,PN}
    megsf = Pair{T,ConstParticalNumberFermionBasis{PN}}[]
    for (coeff,basis) in state.megs 
        for pos in combinations(basis.indices, N÷2)
            k = upper_triangle_count(pos...,h.n)
            for id in h.colptr[k]:h.colptr[k+1]-1
                coeff2, new_basis = apply_SFS_on_basis(h.nzval[id],basis)
                iszero(coeff2) || push!(megsf, Pair(coeff*coeff2, new_basis)) 
            end
        end
    end
    return ConstParticalNumberFermionState(sort_and_combine_terms(megsf))
end

function apply_hamiltonian(h::CSCSimpleFermionHamiltonian, state::ConstParticalNumberFermionState)
    return sort_and_combine_terms(apply_hamiltonian(h.quadratic, state) + apply_hamiltonian(h.quartic, state))
end

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
CPNFB(indices::NTuple) = ConstParticalNumberFermionBasis(indices)
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

function bases_dict(bases::Vector{ConstParticalNumberFermionBasis{PN}}) where PN
    return Dict([b => i for (i,b) in enumerate(bases)])
end

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

function CSCFH_under_bases(h::CSCFH{T,N}, bases::Vector{ConstParticalNumberFermionBasis{PN}},dict::Dict{ConstParticalNumberFermionBasis{PN}, Int64}) where {T,N,PN}
    n = length(bases)
    rowindices, colindices, data = Int[],Int[],T[]
    for base in bases
        for pos in combinations(base.indices, N÷2)
            k = upper_triangle_count(pos...,h.n)
            for id in h.colptr[k]:h.colptr[k+1]-1
                coeff, new_base = apply_SFS_on_basis(h.nzval[id],base)
                if !iszero(coeff) 
                    push!(colindices, dict[base])
                    push!(rowindices, dict[new_base])
                    push!(data, coeff)
                end
            end
        end
    end
    @show rowindices, colindices, data
    return sparse(rowindices, colindices, data, n, n);
end

function CSCSimpleFH_under_bases(h::CSCSimpleFermionHamiltonian, bases::Vector{ConstParticalNumberFermionBasis{PN}},dict::Dict{ConstParticalNumberFermionBasis{PN}, Int64}) where PN
    return CSCFH_under_bases(h.quadratic, bases, dict)+CSCFH_under_bases(h.quartic, bases, dict)
end
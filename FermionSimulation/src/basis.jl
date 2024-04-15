abstract type AbstractBasis end
abstract type AbstractFermionBasis <: AbstractBasis end

# https://en.wikipedia.org/wiki/Fock_space
# PN: particle number
struct FockBasis{PN} <: AbstractFermionBasis
    indices::NTuple{PN,Int}
    function FockBasis(indices::NTuple{PN,Int}) where PN
        new{PN}(TupleTools.sort(indices))
    end
end
Base.:(==)(b1::FockBasis{PN}, b2::FockBasis{PN}) where PN = b1.indices == b2.indices

function Base.show(io::IO,::MIME"text/plain", c::FockBasis{PN}) where PN
    print(io, "|$(join([string(ci) for ci in c.indices], ", "))⟩")
end
Base.show(io::IO, c::FockBasis) = Base.show(io, MIME("text/plain"), c)

function all_bases(PN::Int, n::Int)
    FockBasis.(NTuple{PN, Int}.(combinations(1:n, PN)))
end

function bases_dict(bases::Vector{FockBasis{PN}}) where PN
    return Dict(zip(bases, 1:length(bases)))
end

function CSCFH_under_bases(h::CSCFH{T,N}, bases::Vector{FockBasis{PN}},dict::Dict{FockBasis{PN}, Int64}) where {T,N,PN}
    rowindices, colindices, data = Int[],Int[],T[]
    locations = Tuple.(combinations(1:PN, N÷2))
    _CSCFH_under_bases!(h, bases, dict, rowindices, colindices, data, locations)
end

function _CSCFH_under_bases!(h::CSCFH{T,N}, bases::Vector{FockBasis{PN}},dict::Dict{FockBasis{PN}, Int64}, rowindices, colindices, data, locations) where {T,N,PN}
    n = length(bases)
    for base in bases
        for loc in locations
            pos = getindex.(Ref(base.indices), loc)
            k = upper_triangle_count(pos...,h.n)
            for id in h.colptr[k]:h.colptr[k+1]-1
                fs = h.nzval[id]
                flavors = ntuple(i->fs.ops[i].flavor, N)
                sign = parity_judge(T,flavors,base)
                if !iszero(sign)
                    new_base = apply_SFS_on_basis_withoutsign(flavors,base)
                    push!(colindices, dict[base])
                    push!(rowindices, dict[new_base])
                    push!(data, fs.coeff*sign)
                end
            end
        end
    end
    return sparse(rowindices, colindices, data, n, n)
end

function CSCSimpleFH_under_bases(h::CSCSimpleFermionHamiltonian, bases::Vector{FockBasis{PN}},dict::Dict{FockBasis{PN}, Int64}) where PN
    return CSCFH_under_bases(h.quadratic, bases, dict)+CSCFH_under_bases(h.quartic, bases, dict)
end

function parity_judge(::Type{T},flavors::NTuple{N,Int}, basis::FockBasis{PN}) where {T,N,PN}
    inds = basis.indices
    cflavors = flavors[1:(N÷2)]
    aflavors = flavors[(N÷2+1):N]
    cinds = setdiff(inds, aflavors)
    if !(aflavors ⊆ inds) || !isdisjoint(cflavors, cinds)
        return zero(T)
    end
    return iseven((sum([count(<(aflavors[i]),inds) for i in 1:N÷2]) + sum([count(<(cflavors[i]),cinds) for i in 1:N÷2]))) ? one(T) : -one(T)
end

function apply_SFS_on_basis_withoutsign(flavors::NTuple{N,Int}, basis::FockBasis{PN}) where {N,PN}
    new_indices = replace(basis.indices, Dict(zip(flavors[(N÷2+1):N],flavors[1:(N÷2)]))...)
    return FockBasis(new_indices)
end

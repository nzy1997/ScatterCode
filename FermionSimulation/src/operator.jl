abstract type AbstractOperator end

struct SpaceConfig{D}
    names::NTuple{D, Vector{String}}
end
Base.size(sc::SpaceConfig) = length.(sc.names)
function get_name_of_flavor(sc::SpaceConfig{D}, flavor::Int) where D
    ci = CartesianIndices(size(sc))[flavor]
    return ntuple(k->sc.names[k][ci.I[k]], D)
end

struct FAnnilation{D} <: AbstractOperator
    config::SpaceConfig{D}
    flavor::Int
end
function Base.show(io::IO,::MIME"text/plain", c::FAnnilation)
    print(io, "c($(join(get_name_of_flavor(c.config, c.flavor), ", ")))")
end
Base.show(io::IO, c::FAnnilation) = Base.show(io, MIME("text/plain"), c)

struct Daggered{OT}
    op::OT
end
Base.show(io::IO,::MIME"text/plain", d::Daggered) = print(io, d.op, "†")
Base.show(io::IO, d::Daggered) = Base.show(io, MIME("text/plain"), d)
Base.adjoint(op::Daggered) = op.op
Base.adjoint(op::AbstractOperator) = Daggered(op)
const FCreation{D} = Daggered{FAnnilation{D}}
function FCreation(config::SpaceConfig{D}, flavor::Int) where D
    FAnnilation(config, flavor)'
end

const FermionicOperator = Union{FAnnilation, FCreation}
isfermionic(op::FermionicOperator) = true

# creation operators are moved to the left
struct StandardFermionicString{N, OPT<:NTuple{N, FermionicOperator}} <: AbstractOperator
    ops::OPT
    function StandardFermionicString(ops::OPT) where {N, OPT <: NTuple{N, FermionicOperator}}
        @assert count(x->x isa FCreation, ops) == N ÷ 2 "The number of creation operators should be same as the number of anihilation operators"
        # more sanity checks
        @assert all(i-> (i <= N÷2) ? ops[i] isa FCreation : ops[i] isa FAnnilation, 1:N)
        @assert issorted(ops[1:N÷2], by=x->x.flavor)
        @assert issorted(ops[N÷2+1:N], by=x->x.flavor)
        new{N, OPT}(ops)
    end
end
StandardFermionicString(op1::FermionicOperator, ops::FermionicOperator...) = StandardFermionicString((op1, ops...))

const FSWithCoeff{T} = Pair{T, <:StandardFermionicString}
struct FermionHamiltonian{T<:Number}
    Fs::Vector{FSWithCoeff{T}}
end
function FermionHamiltonian(Fs::Vector{<:Pair{T, <:NTuple{N, FermionicOperator}}}) where {T, N}
    FermionHamiltonian(FSWithCoeff{T}[sort_fops(f.first, f.second) for f in Fs])
end
function sort_fops(coeff, ops::NTuple{N, FermionicOperator}) where N
    perm = TupleTools.sortperm(ops, lt=(x, y)->(x isa FCreation && (y isa FAnnilation || x.flavor < y.flavor)) || (x isa FAnnilation && y isa FAnnilation && x.flavor < y.flavor))
    coeff = isevenperm(perm) ? coeff : -coeff
    return coeff => StandardFermionicString(getindex.(Ref(ops), perm))
end

function Base.push!(h::FermionHamiltonian, p::Pair{T, <:Tuple}) where {T}
    item = sort_fops(p.first, p.second)
    push!(h.Fs, item)
    return h
end

# credit: https://github.com/toivoh/ConicHulls.jl/blob/888109ab6fe7d334ffa48e649947f0c293a967ae/src/Common.jl#L37
function isevenperm(p)
    @assert isperm(p)
    n = length(p)
    used = falses(n)
    even = true
    for k = 1:n
        if used[k]; continue; end
        # Each even cycle flips even (an odd number of times)
        used[k] = true
        j = p[k]
        while !used[j]
            used[j] = true
            j = p[j]
            even = !even
        end
    end
    even
end

Base.:(+)(h1::FermionHamiltonian{T}, h2::FermionHamiltonian{T}) where T = FermionHamiltonian(vcat(h1.Fs, h2.Fs))
# Base.:(*)(j0::T, fs::StandardFermionicString) where T = FermionHamiltonian(Pair{T, <:StandardFermionicString}[j0 => fs])
Base.:(-)(h1::FermionHamiltonian{T}) where T = FermionHamiltonian(FSWithCoeff{T}[-x.first=>x.second for x in h1.Fs])
Base.:(-)(h1::FermionHamiltonian{T}, h2::FermionHamiltonian{T}) where T = h1 + (-h2)
# Base.:(-)(fs::StandardFermionicString) = - 1.0 * fs


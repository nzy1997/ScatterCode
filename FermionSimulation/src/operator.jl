abstract type AbstractOperator end

struct Fop <: AbstractOperator
	iscreation::Bool
	flavor::Int
end
function Base.show(io::IO, ::MIME"text/plain", c::Fop)
	print(io, "c$(c.flavor)$(c.iscreation ? "†" : "")")
end
Base.show(io::IO, c::Fop) = Base.show(io, MIME("text/plain"), c)
creation(flavor::Int) = Fop(true, flavor)
annilation(flavor::Int) = Fop(false, flavor)

# creation operators are moved to the left
struct StandardFermionicString{T, N} <: AbstractOperator
	coeff::T
	ops::NTuple{N, Fop}
	function StandardFermionicString(coeff::T, ops::NTuple{N, Fop}) where {T, N}
		@assert count(x -> x.iscreation, ops) == N ÷ 2 "The number of creation operators should be same as the number of anihilation operators"
		coeff_, ops_ = sort_fops(coeff, ops)
		new{T, N}(coeff_, ops_)
	end
end
function sort_fops(coeff, ops::NTuple{N, Fop}) where N
	perm = TupleTools.sortperm(ops, lt = (x, y) -> (x.iscreation && (!y.iscreation || x.flavor < y.flavor)) || (!x.iscreation && !y.iscreation && x.flavor < y.flavor))
	coeff = isevenperm(perm) ? coeff : -coeff
	return coeff, getindex.(Ref(ops), perm)
end
StandardFermionicString(coeff::Number, op1::Fop, ops::Fop...) = StandardFermionicString(coeff, (op1, ops...))
Base.:*(x::Number, fs::StandardFermionicString) = StandardFermionicString(x * fs.coeff, fs.ops)
Base.:-(fs::StandardFermionicString) = (-1) * fs

struct SimpleFermionHamiltonian{T <: Number}
	quadratic::Vector{StandardFermionicString{T, 2}}
	quartic::Vector{StandardFermionicString{T, 4}}
	function SimpleFermionHamiltonian(quadratic::Vector{StandardFermionicString{T, 2}}) where T
		SimpleFermionHamiltonian(quadratic, StandardFermionicString{T, 4}[])
	end
	function SimpleFermionHamiltonian(quadratic::Vector{StandardFermionicString{T, 4}}) where T
		SimpleFermionHamiltonian(StandardFermionicString{T, 2}[], quadratic)
	end
	function SimpleFermionHamiltonian(quadratic::Vector{StandardFermionicString{T, 2}}, quartic::Vector{StandardFermionicString{T, 4}}) where T
		new{T}(quadratic, quartic)
	end
end
function Base.push!(h::SimpleFermionHamiltonian, p::StandardFermionicString{T, 2}) where {T}
	push!(h.quadratic, p)
	return h
end
function Base.push!(h::SimpleFermionHamiltonian, p::StandardFermionicString{T, 4}) where {T}
	push!(h.quartic, p)
	return h
end

Base.:(+)(h1::SimpleFermionHamiltonian{T}, h2::SimpleFermionHamiltonian{T}) where T = SimpleFermionHamiltonian(vcat(h1.quadratic, h2.quadratic), vcat(h1.quartic, h2.quartic))
Base.:*(x::Number, h::SimpleFermionHamiltonian{T}) where T = SimpleFermionHamiltonian(x .* h.quadratic, x .* h.quartic)
Base.:(-)(h1::SimpleFermionHamiltonian{T}) where T = (-1) * h1
Base.:(-)(h1::SimpleFermionHamiltonian{T}, h2::SimpleFermionHamiltonian{T}) where T = h1 + (-h2)

matrix_size(N::Int,n::Int) = binomial(n,N ÷ 2)

struct CSCFH{T,N}
	n::Int
	colptr::Vector{Int}
	rowval::Vector{Int}
    nzval::Vector{StandardFermionicString{T, N}}
	function CSCFH{T, N}(n::Int, colptr::Vector{Int}, rowval::Vector{Int}, nzval::Vector{StandardFermionicString{T, N}}) where {T,N}
		@assert length(colptr) == matrix_size(N,n) + 1
		@assert length(rowval) == length(nzval) == colptr[end] - 1
		new{T,N}(n, colptr, rowval, nzval)
	end
end


function get_ac_index(fs::StandardFermionicString{T, 2},n::Int) where T
    return fs.ops[1].flavor, fs.ops[2].flavor
end

function get_ac_index(fs::StandardFermionicString{T, 4},n::Int) where T
    k1 = upper_triangle_count(fs.ops[1].flavor, fs.ops[2].flavor, n)
    k2 = upper_triangle_count(fs.ops[3].flavor, fs.ops[4].flavor, n)
    return k1, k2
end

function sort_SFS(fs::StandardFermionicString{T, N}, n::Int) where {T,N}
    k1,k2 = get_ac_index(fs,n)
    return k1 + matrix_size(N,n) * (k2 - 1)
end

function make_sfs(k1::Int,k2::Int,N::Int,coeff::T,n::Int) where {T}
    if N == 2
        return StandardFermionicString(coeff, (creation(k1), annilation(k2)))
    elseif N == 4
        kc1,kc2=upper_triangle_count_inv(k1, n)
        ka1,ka2 = upper_triangle_count_inv(k2, n)
        return StandardFermionicString(coeff, (creation(kc1), creation(kc2), annilation(ka1), annilation(ka2)))
    end
end

function CSCFH(quadratic::Vector{StandardFermionicString{T, N}}, n::Int) where {T,N}
	# sort the COO matrix by column
    ms = matrix_size(N,n)
	order = sortperm(1:length(quadratic); by = i -> sort_SFS(quadratic[i], n))
	colptr, rowval = zeros(Int, ms + 1), zeros(Int, length(quadratic))
    nzval = StandardFermionicString{T, N}[]
	k = 0
	ipre, jpre = 0, 0
	colptr[1] = 1
	for idx in order
        i, j = get_ac_index(quadratic[idx],n)
		v = quadratic[idx].coeff
		# values with the same indices are accumulated
		if i == ipre && j == jpre
			v += nzval[k].coeff
            pop!(nzval)
            push!(nzval, make_sfs(i,j,N,v,n))
		else
			k += 1
			if j != jpre
				# a new column starts
				colptr[jpre+1:j+1] .= k
			end
			rowval[k] = i
			push!(nzval, make_sfs(i,j,N,v,n))
			ipre, jpre = i, j
		end
	end
	colptr[jpre+1:end] .= k + 1
	resize!(rowval, k)
	return CSCFH{T,N}(n, colptr, rowval, nzval)
end

struct CSCSimpleFermionHamiltonian{T}
	quadratic::CSCFH{T,2}
	quartic::CSCFH{T,4}
end

function CSCSimpleFermionHamiltonian(h::SimpleFermionHamiltonian{T},n::Int) where T
    CSCSimpleFermionHamiltonian(CSCFH(h.quadratic, n), CSCFH(h.quartic, n))
end

# struct SpaceConfig{D}
#     names::NTuple{D, Vector{String}}
# end
# Base.size(sc::SpaceConfig) = length.(sc.names)
# function get_name_of_flavor(sc::SpaceConfig{D}, flavor::Int) where D
#     ci = CartesianIndices(size(sc))[flavor]
#     return ntuple(k->sc.names[k][ci.I[k]], D)
# end

# function print_fop_with_sc(io::IO, ::FAnnilation, config::SpaceConfig)
#     print(io, "c($(join(get_name_of_flavor(config, config), ", ")))")
# end

# #### No need
# struct Daggered{OT}
#     op::OT
# end
# Base.show(io::IO,::MIME"text/plain", d::Daggered) = print(io, d.op, "†")
# Base.show(io::IO, d::Daggered) = Base.show(io, MIME("text/plain"), d)
# Base.adjoint(op::Daggered) = op.op
# Base.adjoint(op::AbstractOperator) = Daggered(op)
# const FCreation = Daggered{FAnnilation}
# function FCreation(flavor::Int)
#     FAnnilation(flavor)'
# end

# const FermionicOperator = Union{FAnnilation, FCreation}
# isfermionic(op::FermionicOperator) = true


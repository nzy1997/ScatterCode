# required interfaces: nsite
abstract type AbstractModel end

# XY model
struct XYModel <: AbstractModel
    Jxy::Float64   # coupling strength
    h::Float64  # magnetic field strength
    n::Int  # number of sites
end
# get the number of sites
nsite(m::XYModel) = m.n

# Fermionic quadratic Hamiltonian for solving XY model exactly
function fermionic_quadratic(model::XYModel)
    n = nsite(model)
    Jxy, h = model.Jxy, model.h
    A=zeros(n, n)
    A[2:(n+1):end] .= -Jxy/2
    A[(n+1):(n+1):end] .= -Jxy/2
    A[1:(n+1):end] .= h
    return A
end

# ```math
# H = H_0 + H_1 \cos(\omega t)
# ```
struct DrivedHamiltonian{B1<:AbstractBlock, B2<:AbstractBlock}
    H0::B1  # static Hamiltonian
    H1::B2  # driving Hamiltonian
    frequency::Float64  # frequency of the driving term
end
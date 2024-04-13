abstract type AbstractBasis end
abstract type AbstractFermionBasis <: AbstractBasis end

struct ConstParticalNumberFermionBasis{PN} <: AbstractFermionBasis
    sc::SpaceConfig
    indices::NTuple{PN,Int}
end

function Base.show(io::IO,::MIME"text/plain", c::ConstParticalNumberFermionBasis{PN}) where PN
    [print(io, "$(get_name_of_flavor(c.sc,c.indices[i]))") for i in 1:PN]
end
Base.show(io::IO, c::ConstParticalNumberFermionBasis) = Base.show(io, MIME("text/plain"), c)

function all_bases(sc::SpaceConfig, PN::Int)
    bases = ConstParticalNumberFermionBasis{PN}[]
    for pos in combinations(1:reduce(*,size(sc)), PN)
        push!(bases, ConstParticalNumberFermionBasis{PN}(sc,(pos...,)))
    end
    return bases
end


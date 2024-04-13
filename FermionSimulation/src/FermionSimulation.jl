module FermionSimulation
using TupleTools
using Combinatorics
using Graphs
using LinearAlgebra

# basis
export ConstParticalNumberFermionBasis, all_bases

# operator
export SpaceConfig, FermionOperator, FAnnilation, FCreation, isfermionic, StandardFermionicString,get_name_of_flavor
export FermionHamiltonian

include("operator.jl")
include("basis.jl")

end

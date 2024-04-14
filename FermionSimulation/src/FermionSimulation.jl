module FermionSimulation
using TupleTools
using Combinatorics
using Graphs
using LinearAlgebra

# operator
export Fop, StandardFermionicString, creation, annilation
export SimpleFermionHamiltonian, CSCFH, CSCSimpleFermionHamiltonian
export upper_triangle_count,upper_triangle_count_inv

# basis
export ConstParticalNumberFermionBasis, all_bases, ConstParticalNumberFermionState ,sort_and_combine_terms
export apply_SFS_on_basis, apply_hamiltonian

include("utils.jl")
include("operator.jl")
include("basis.jl")

end

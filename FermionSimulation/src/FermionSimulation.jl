module FermionSimulation
using TupleTools
using Combinatorics
using Graphs
using LinearAlgebra
using SparseArrays

# operator
export Fop, StandardFermionicString, creation, annilation
export SimpleFermionHamiltonian, CSCFH, CSCSimpleFermionHamiltonian
export upper_triangle_count,upper_triangle_count_inv

# basis
export ConstParticalNumberFermionBasis,CPNFB, all_bases,bases_dict, ConstParticalNumberFermionState ,sort_and_combine_terms
export apply_SFS_on_basis, CSCFH_under_bases, CSCSimpleFH_under_bases

include("utils.jl")
include("operator.jl")
include("basis.jl")

end

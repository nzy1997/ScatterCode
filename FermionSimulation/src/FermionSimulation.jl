module FermionSimulation
using TupleTools
using Combinatorics
using Graphs
using LinearAlgebra
using SparseArrays
using KrylovKit
using Scattering
using Scattering:GraphWithTails, gaussian_packet_ontail

# operator
export Fop, StandardFermionicString, creation, annilation
export SimpleFermionHamiltonian, CSCFH, CSCSimpleFermionHamiltonian
export upper_triangle_count,upper_triangle_count_inv

# basis
export FockBasis, all_bases,bases_dict, ConstParticalNumberFermionState ,sort_and_combine_terms
export apply_SFS_on_basis, CSCFH_under_bases, CSCSimpleFH_under_bases
export kron_state, identical_particle_wave

# structure
export get_hamiltonian_from_sg, get_CSCSimpleFH, example_gt_cz, example_gt_h, example_gt_doublechain

# time_evo
export time_evolve, expect_value, position_operators, simulate_gt

include("utils.jl")
include("operator.jl")
include("basis.jl")
include("structure.jl")
include("time_evolve.jl")
end

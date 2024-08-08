module Scattering
using Graphs
using Optim
using LineSearches
using ForwardDiff
using LinearAlgebra

using FFTW
using Yao
using LuxorGraphPlot
using LuxorGraphPlot.Luxor
using SparseArrays

# scatter_center
export ScatterGraph, abs_loss, perm_adjacency, matrix_adjacency
export scatter_matrix, get_u1, get_u2, get_r1,relection_rate

# structures
export andrew_phase_shift, andrew_basis_change, andrew_momentum_filter, andrew_momentum_separator, andrew_hadamard, wuzi

# optimization
export OptimizedResult, train_loss, optimize_momentum, optimize_weighted_momentum, weighted_gra!

# universal_check
export fraction_approximate, unitary_decomposition, universal_check,random_unitary

# simulate
export graph_with_tails, center_graph, weight_oncenter, weight_ontail
export get_hamiltonian
export generate_gaussian_packet
export simulate_chain, simulate_graph_with_tails, simulate_ScatterGraph

# visualize
export animate_wave, plot_transmission, animate_graph

include("scattering_center.jl")
include("structures.jl")
include("optimization.jl")
include("universal_check.jl")

include("simulate.jl")
include("visualize.jl")
end

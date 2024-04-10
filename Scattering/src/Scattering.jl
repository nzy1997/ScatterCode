module Scattering
using Graphs
using Optim
using ForwardDiff
using LinearAlgebra

using FFTW
using Yao
using LuxorGraphPlot
using LuxorGraphPlot.Luxor
using SparseArrays

# scatter_center
export ScatterGraph, abs2_loss, perm_adjacency, scatter_matrix, matrix_adjacency

# optimization
export train_loss, optimize_momentum, optimize_weighted_momentum, weighted_gra!

# universal_check
export fraction_approximate, unitary_decomposition, universal_check,random_unitary

# simulate
export graph_with_tails, center_graph, weight_oncenter, weight_ontail
export get_hamiltonian
export generate_gaussian_packet
export simulate_chain, simulate_graph_with_tails

# visualize
export animate_wave, plot_transmission

include("scattering_center.jl")
include("optimization.jl")
include("universal_check.jl")

include("simulate.jl")
include("visualize.jl")
end

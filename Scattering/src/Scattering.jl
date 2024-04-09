module Scattering
using Graphs
using SimpleWeightedGraphs
using Enzyme
using Optim
using ForwardDiff


# scatter_center
export ScatterGraph, abs2_loss, perm_adjacency, scatter_matrix, matrix_adjacency

# optimization
export train_loss, optimize_momentum, optimize_weighted_momentum, weighted_gra!

# universal_check
export fraction_approximate, unitary_decomposition, universal_check,random_unitary

include("scattering_center.jl")
include("optimization.jl")
include("universal_check.jl")
end

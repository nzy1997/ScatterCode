using FermionSimulation
using Test

@testset "operator.jl" begin
    include("operator.jl")
end

@testset "basis.jl" begin
    include("basis.jl")
end

@testset "graph_struct.jl" begin
    include("graph_struct.jl")
end

@testset "time_evolve.jl" begin
    include("time_evolve.jl")
end
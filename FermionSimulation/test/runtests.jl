using FermionSimulation
using Test

@testset "operator.jl" begin
    include("operator.jl")
end

@testset "basis.jl" begin
    include("basis.jl")
end

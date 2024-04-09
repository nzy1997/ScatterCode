using Scattering
using Test

@testset "scatter_center" begin
    include("scatter_center.jl")
end

@testset "optimization" begin
    include("optimization.jl")
end

@testset "universal_check" begin
    include("universal_check.jl")
end
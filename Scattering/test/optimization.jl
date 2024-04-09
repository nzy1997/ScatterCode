using Test
using Scattering
using Scattering.Graphs
using LinearAlgebra
import Random
Random.seed!(1234)

@testset "gradient" begin
    g = ScatterGraph(SimpleGraph(Edge.([(1,5), (2,7), (3, 8), (4, 10), (5,6), (5,8), (6,7), (7,10), (8,9), (9, 10)])), [1, 2], [3, 4])
    G = zeros(11)
    x0 = randn(11)
    δ = 1e-8
    direction = rand(11)
    direction = direction/norm(direction)
    x1 = x0 + δ * direction
    gra = weighted_gra!(G, x0, g)
    @test abs2_loss(g, exp(pi * im * x1[1]), x1[2:end]) ≈ abs2_loss(g, exp(pi * im * x0[1]), x0[2:end]) +  δ * dot(gra, direction)
end

@testset "optimize_momentum" begin
    g = ScatterGraph(SimpleGraph(Edge.([(1,5), (2,7), (3, 8), (4, 10), (5,6), (5,8), (6,7), (7,10), (8,9), (9, 10)])), [1, 2], [3, 4])
    x, min = optimize_momentum(g)
    @test isapprox(min, 0.0, atol=1e-20)
    @test isapprox(x[1], 0.25, atol=1e-20)
end

@testset "optimize_weighted_momentum" begin
    g = ScatterGraph(SimpleGraph(Edge.([(1,5), (2,7), (3, 8), (4, 10), (5,6), (5,8), (6,7), (7,10), (8,9), (9, 10)])), [1, 2], [3, 4])
    x0 = rand(Float64, length(edges(g.graph))+1)
    x, min = optimize_weighted_momentum(g,x0)
    @show x, min
    @test isapprox(min, 0.0, atol=1e-6)
end

@testset "optimize_weighted_momentum wuzi" begin
    g = ScatterGraph(SimpleGraph(Edge.([(1,5), (5,7), (3, 5), (4, 9), (7,8), (4,8), (6,7), (7,9), (2,9)])), [1, 2], [3, 4])
    x0 = rand(Float64, length(edges(g.graph))+1)
    x, min = optimize_weighted_momentum(g, x0)
    z=exp(pi*im*x[1])
    s=scatter_matrix(g, z, x[2:end])
    @test isapprox(s[1:2,1:2], zeros(2,2), atol=1e-5)
    @test isapprox(s[3:4,3:4], zeros(2,2), atol=1e-5)
end

@testset "optimize_weighted_momentum hardmard" begin
    g = ScatterGraph(SimpleGraph(Edge.([(1,7), (5,6),(4,6),(3,11),(6,7),(6,10),(5,9),(5,8),(7,8),(7,10),(8,9),(9,10),(8,11),(9,13),(10,11),(10,13),(8,12),(11,12),(2,12),(12,13)])), [1, 2], [3, 4])
    x0 = rand(Float64, length(edges(g.graph))+1)
    x, min = optimize_weighted_momentum(g, x0)
    @show x,min
    z=exp(pi*im*x[1])
    s=scatter_matrix(g, z, x[2:end])
    @test isapprox(s[1:2,1:2], zeros(2,2), atol=1e-5)
    @test isapprox(s[3:4,3:4], zeros(2,2), atol=1e-5)
    u1=s[1:2,3:4]
    u2=s[3:4,1:2]
    ulist = [u1, u2]
    @show universal_check(ulist,1)
end
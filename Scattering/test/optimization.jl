using Test
using Scattering
using Scattering.Graphs
using LinearAlgebra
import Random

@testset "gradient" begin
    g = andrew_basis_change()
    G = zeros(11)
    x0 = randn(11)
    δ = 1e-8
    direction = rand(11)
    direction = direction/norm(direction)
    x1 = x0 + δ * direction
    gra = weighted_gra!(G, x0, g)
    @test abs_loss(g.graph, x1[2:end],g.input_vertex,g.output_vertex, exp(pi * im * x1[1])) ≈ abs_loss(g.graph, x0[2:end],g.input_vertex,g.output_vertex, exp(pi * im * x0[1])) +  δ * dot(gra, direction)
end

@testset "optimize_momentum" begin
    g = andrew_basis_change()
    Random.seed!(123523)
    x, min = optimize_momentum(g,rand())
    @test isapprox(min, 0.0, atol=1e-10)
    @test isapprox(x[1], 0.25, atol=1e-8)
end

@testset "optimize_weighted_momentum" begin
    Random.seed!(123523)
    x0 = rand(Float64, 10)
    g = andrew_basis_change(x0)
    outcome = optimize_weighted_momentum(g,rand(Float64))
    @test isapprox(outcome.minimum, 0.0, atol=1e-6)
end

@testset "optimize_weighted_momentum wuzi" begin
    Random.seed!(123)
    x0 = rand(Float64, 9)
    g = wuzi(x0)
    outcome= optimize_weighted_momentum(g,rand(Float64))
    s = outcome.s
    @test isapprox(s[1:2,1:2], zeros(2,2), atol=1e-5)
    @test isapprox(s[3:4,3:4], zeros(2,2), atol=1e-5)
end

@testset "optimize_weighted_momentum hardmard" begin
    Random.seed!(435)
    x0 = rand(Float64, 20)
    g = andrew_hadamard(x0)
    outcome = optimize_weighted_momentum(g,rand(Float64))
    s=outcome.s
    @test isapprox(s[1:2,1:2], zeros(2,2), atol=1e-5)
    @test isapprox(s[3:4,3:4], zeros(2,2), atol=1e-5)
    u1=get_u1(s)
    u2=get_u2(s)
    ulist = [u1, u2]
    @test universal_check(ulist,1)
end
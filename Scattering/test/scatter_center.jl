using Scattering
using Test
using Scattering.Graphs

@testset "abs2_loss" begin
    g = ScatterGraph(SimpleGraph(Edge.([(1,5), (2,7), (3, 8), (4, 10), (5,6), (5,8), (6,7), (7,10), (8,9), (9, 10)])), [1, 2], [3, 4])
    z = exp(-pi/4*im)
    @test isapprox(abs2_loss(g, z),0.0, atol=1e-20)
end

@testset "perm_adjacency" begin
    g = ScatterGraph(SimpleGraph(Edge.([(1,5), (2,7), (3, 8), (4, 10), (5,6), (5,8), (6,7), (7,10), (8,9), (9, 10)])), [1, 5], [6, 4])
    W = matrix_adjacency(g.graph, ones(Float64, ne(g.graph)))
    @test perm_adjacency(g,W)[1,2] == 1
    @test perm_adjacency(g,W)[3,7] == 1
end

@testset "scatter_matrix" begin
    g = ScatterGraph(SimpleGraph(Edge.([(1,5), (2,7), (3, 8), (4, 10), (5,6), (5,8), (6,7), (7,10), (8,9), (9, 10)])), [1, 2], [3, 4])
    z = exp(-pi/4*im)
    s=scatter_matrix(g, z)
    @test isapprox(s[1:2,1:2], zeros(2,2), atol=1e-10)
    @test isapprox(s[3:4,3:4], zeros(2,2), atol=1e-10)
    U=-1/sqrt(2)*[im 1; 1 im]
    @test isapprox(s[1:2,3:4], U, atol=1e-10)
    @test isapprox(s[3:4,1:2], transpose(U), atol=1e-10)
end

@testset "scatter_center hardmard" begin
    g = ScatterGraph(SimpleGraph(Edge.([(1,7), (5,6),(4,6),(3,11),(6,7),(6,10),(5,9),(5,8),(7,8),(7,10),(8,9),(9,10),(8,11),(9,13),(10,11),(10,13),(8,12),(11,12),(2,12),(12,13)])), [1, 2], [3, 4])
    z = exp(-pi/2*im)
    @show abs2_loss(g, z)
    @test isapprox(abs2_loss(g, z),0.0, atol=1e-20)
end
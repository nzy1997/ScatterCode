using Test, Scattering, Graphs, Yao

@testset "chain" begin
    n = 1000
    x0 = 100
    k0 = π/4
    Δx = 20
    # gaussian packet
    wave = generate_gaussian_packet(n, x0; k0, Δx)
    m, s = Scattering.meanstd(abs.(wave))
    @test isapprox(m, x0; atol=1e-1)
    @test isapprox(s, Δx; atol=1e-1)

    # chain
    α = 1.0
    v = simulate_chain(n; x0, k0, Δx, Δt=1.0, Nt=1000, α)
    p2 = sum(abs2, v[n÷2+1:n])
    T2 = 1/(1+α^2 / sin(k0)^2)
    println("T² = $(1/(1+α^2 / sin(k0)^2)), got $(p2)")
    @test isapprox(T2, p2; atol=1e-2)
end

@testset "graph with tails" begin
    α = 1.0
    gt = graph_with_tails(path_graph(3); heads=[1, 3], tailsize=499)
    h0 = Diagonal([0, α, 0])
    cache = Vector{ComplexF64}[]
    v = simulate_graph_with_tails(gt, h0; Nt=1000, cache)
    @test isapprox(weight_ontail(v, gt, 1), 2/3; atol=1e-2)
    @test isapprox(weight_ontail(v, gt, 2), 1/3; atol=1e-2)
    @test isapprox(weight_oncenter(v, gt), 0; atol=1e-2)
    @test length(cache) == 1000
end

@testset "graph with tails 2" begin
    gt = graph_with_tails(SimpleGraph(Edge.([(1,3), (2,3), (3,4), (4,5), (5,6), (6,7), (7,4)])); heads=[1, 2], tailsize=499)
    h0 = Diagonal(zeros(7))
    v = simulate_graph_with_tails(gt, h0; Nt=1000)
    @test isapprox(weight_ontail(v, gt, 1), 0; atol=1e-2)
    @test isapprox(weight_ontail(v, gt, 2), 1.0; atol=1e-2)
    @test isapprox(weight_oncenter(v, gt), 0; atol=1e-2)
end
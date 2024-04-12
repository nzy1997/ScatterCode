using Test, Scattering, Graphs, LinearAlgebra
using LuxorGraphPlot: Luxor

@testset "animate wave" begin
    α = 1.0
    gt = graph_with_tails(SimpleGraph(Edge.([(1,3), (2,3), (3,4), (4,5), (5,6), (6,7), (7,4)])); heads=[1, 2], tailsize=499)
    h0 = Diagonal(zeros(7))
    waves = Vector{ComplexF64}[]
    v = simulate_graph_with_tails(gt, h0; Nt=1000, cache=waves)
    @test animate_wave(gt, waves, step=100) isa Luxor.AnimatedGif
    # @test plot_transmission(andrew_momentum_separator(); outtail=2) isa Luxor.Drawing
end
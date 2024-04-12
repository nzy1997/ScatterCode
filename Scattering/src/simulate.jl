struct GraphWithTails
    graph::SimpleGraph{Int}
    center::Vector{Int}
    tails::Vector{Vector{Int}}
end
graph_with_tails(g::ScatterGraph; heads::Vector{Int}, tailsize::Int) = graph_with_tails(g.graph; heads, tailsize)

function graph_with_tails(g0::SimpleGraph; heads::Vector{Int}, tailsize::Int)
    ntail = length(heads)
    nv0 = nv(g0)
    graph = SimpleGraph(tailsize * ntail + nv(g0))
    for edge in edges(g0)
        add_edge!(graph, edge)
    end
    tails = Vector{Int}[]
    for i=1:ntail
        tl = [heads[i], [nv0+(i-1)*tailsize+j for j=1:tailsize]...]
        push!(tails, tl)
        for j=1:tailsize
            add_edge!(graph, tl[j], tl[j+1])
        end
    end
    return GraphWithTails(graph, collect(1:nv0), tails)
end

center_graph(gt::GraphWithTails) = induced_subgraph(gt.graph, gt.center)[1]
weight_ontail(wave::AbstractVector, gt::GraphWithTails, index::Int) = sum(v->abs2(wave[v]), gt.tails[index])
weight_oncenter(wave::AbstractVector, gt::GraphWithTails) = sum(v->abs2(wave[v]), gt.center)

function get_hamiltonian(gt::GraphWithTails, h0::AbstractMatrix, ::Type{T}=ComplexF64) where T
    h = get_hamiltonian(gt.graph, T)
    h[gt.center, gt.center] += h0
    return h
end

function get_hamiltonian(graph::SimpleGraph, ::Type{T}=ComplexF64) where T
    return one(T)/2 * adjacency_matrix(graph)
end

function generate_gaussian_packet(n::Int, x0; k0, Δx)
    @assert -π <= k0 <= π "momentum should be in range [-π, π], got $k0"
    Δk = 1/Δx
    ks = 2π .* (collect(0:n-1) .- (n÷2)) ./ n
    packet_k_half = exp.((ks .- k0) .^ 2 ./ (-2*Δk^2))
    packet_k =[packet_k_half[n÷2+1:end]... ,packet_k_half[1:n÷2]...]
    amplitude = circshift!(ifft(packet_k), x0)
    return normalize!(amplitude)
end

function meanstd(probs, xs = 0:length(probs)-1)
    probs = normalize!(probs, 1)
    n = length(probs)
    mean = sum(probs .* xs)
    return mean, sqrt(sum(i-> probs[i]*(xs[i] - mean)^2, 1:n))
end

function simulate_chain(n::Int;
        x0::Int=round(Int, 0.1*n), Δx::Real=0.02*n,
        k0::Real=π/4, Δt::Real=1.0, Nt::Int=1000, α::Real=0.0,
        cache::Union{Vector, Nothing}=nothing
        )
    wave = generate_gaussian_packet(n, x0; k0, Δx)
    h = get_hamiltonian(path_graph(n))
    h[n÷2, n÷2] += α
    return simulate!(wave, h; Nt, Δt, cache)
end

function simulate!(wave, hmat;
    Nt, Δt, cache::Union{Vector, Nothing}=nothing)
    n = length(wave)
    reg = arrayreg(wave; nlevel=n)
    block = matblock(hmat; nlevel=size(hmat, 1))
    te = time_evolve(block, Δt)
    for _ = 1:Nt
        apply!(reg, te)
        cache === nothing || push!(cache, copy(statevec(reg)))
    end
    return statevec(reg)
end

function gaussian_packet_ontail(gt::GraphWithTails; intail::Int, k0, x0, Δx)
    packet = generate_gaussian_packet(length(gt.tails[intail]), x0; k0, Δx)
    v = zeros(ComplexF64, nv(gt.graph))
    v[gt.tails[intail]] .= packet # [end:-1:1]
    return v
end

function simulate_graph_with_tails(gt::GraphWithTails, h0::AbstractMatrix;
        intail::Int=1, k0::Real=π/4,
        x0::Int=round(Int, 0.8*length(gt.tails[intail])),
        Δx::Real=0.04*length(gt.tails[intail]), 
        Nt::Int=1000, Δt::Real=1.0,
        cache::Union{Vector, Nothing}=nothing
        )
    wave = gaussian_packet_ontail(gt; intail, k0, x0, Δx)
    h = get_hamiltonian(gt, h0)
    simulate!(wave, h; Nt, Δt, cache)
end

function simulate_ScatterGraph(g::ScatterGraph, x1::Real ;Nt=1500, tailsize::Int=500, intail::Int=1) # x1 is k/pi where k is momentum
    gt = graph_with_tails(g; heads=[g.input_vertex...,g.output_vertex...], tailsize)
    h0 = matrix_adjacency(g.graph, g.weights .- 1) ./ 2  # difference in weight
    k0 = (mod(x1,2.0) > 1.0 ? (mod(x1,2.0)-2.0)*pi : mod(x1,2.0)*pi)
    waves = Vector{ComplexF64}[]
    simulate_graph_with_tails(gt, h0;k0, Nt, cache=waves,intail)
    return waves,gt
end
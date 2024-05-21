struct TimeEvolver{T, MT <: AbstractMatrix{T}, VT <: AbstractVector{T}}
    H::MT
    psi::VT
end
function evolve(evolver::TimeEvolver{T, MT}, tstep::Float64) where {T, MT}
    psi, info = exponentiate(evolver.H, -im * tstep, evolver.psi,
        krylovdim = min(1000, size(evolver.H, 1)),
        ishermitian = true,
        eager = true)
    return TimeEvolver(evolver.H, psi)
end

function time_evolve(psi::AbstractVector{T}, H::AbstractMatrix{T}, Nt::Int, tstep::Float64; cache::Union{Vector, Nothing} = nothing, operators::AbstractVector = []) where {T}
    cache === nothing || push!(cache, psi)
    expects = [[expect_value(psi, operator) for operator in operators]]
    evo = TimeEvolver(H, psi)
    for i ∈ 1:Nt
        evo = evolve(evo, tstep)
        cache === nothing || push!(cache, evo.psi)
        @show i
        push!(expects, [expect_value(evo.psi, operator) for operator in operators])
    end
    return expects
end

function expect_value(psi::AbstractVector{T}, h::AbstractMatrix{T}) where {T}
    return real(psi' * h * psi)
end

function position_operators(n::Int, bases::Vector{FockBasis{PN}}, dict::Dict{FockBasis{PN}, Int64}) where PN
    return [ComplexF64.(CSCSimpleFH_under_bases(CSCSimpleFermionHamiltonian(SimpleFermionHamiltonian([StandardFermionicString(1.0, creation(i), annilation(i))]), n), bases, dict)) for i in 1:n]
end

function simulate_gt(gt::GraphWithTails, k0::Float64, intails::Vector{Int};Δt::Real=1.0, Nt::Int=999,tag::Bool=true,hs::Union{Nothing,SimpleFermionHamiltonian}=nothing)
    particle_num = length(intails)
    wave = kron_state([gaussian_packet_ontail(gt; intail, k0, x0 = round(Int, 0.2 * length(gt.tails[intail])), Δx = 0.02 * length(gt.tails[intail])) for intail in intails]...)
    wave=wave/norm(wave)
    H ,bases, dict,n = get_CSCSimpleFH(gt.graph, particle_num)
    if hs !== nothing
        H = H + CSCSimpleFH_under_bases(CSCSimpleFermionHamiltonian(hs, n), bases, dict)
    end
    if tag
        operators = position_operators(n, bases, dict)
        waves = time_evolve(wave, ComplexF64.(H), Nt, Δt; operators)
    else
        waves = Vector{ComplexF64}[]
        time_evolve(wave, ComplexF64.(H), Nt, Δt; cache = waves)
    end
    return waves
end

function simulate_gt(gt::GraphWithTails, k::Vector{Float64}, intails::Vector{Int};Δt::Real=1.0, Nt::Int=999,tag::Bool=true,hs::Union{Nothing,SimpleFermionHamiltonian}=nothing)
    particle_num = length(intails)
    wave = kron_state([gaussian_packet_ontail(gt; intail=intails[i], k0 = k[i], x0 = round(Int, 0.2 * length(gt.tails[intails[i]])), Δx = 0.02 * length(gt.tails[intails[i]])) for i in 1:particle_num]...)
    wave = wave/norm(wave)
    H ,bases, dict,n = get_CSCSimpleFH(gt.graph, particle_num)
    if hs !== nothing
        H = H + CSCSimpleFH_under_bases(CSCSimpleFermionHamiltonian(hs, n), bases, dict)
    end
    if tag
        operators = position_operators(n, bases, dict)
        waves = time_evolve(wave, ComplexF64.(H), Nt, Δt; operators)
    else
        waves = Vector{ComplexF64}[]
        time_evolve(wave, ComplexF64.(H), Nt, Δt; cache = waves)
    end
    return waves
end
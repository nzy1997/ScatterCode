struct TimeEvolver{T, MT<:AbstractMatrix{T}, VT<:AbstractVector{T}}
    H::MT
    psi::VT
end
function evolve(evolver::TimeEvolver{T, MT}, tstep::Float64) where {T, MT}
    psi, info = exponentiate(evolver.H, -im * tstep, evolver.psi,
        krylovdim=min(1000, size(evolver.H, 1)),
        ishermitian=true,
        eager=true)
    return TimeEvolver(evolver.H, psi)
end

function time_evolve(psi::AbstractVector{T}, H::AbstractMatrix{T}, Nt::Int, tstep::Float64) where {T}
    psis = [psi]
    evo = TimeEvolver(H, psi)
    for _ = 1:Nt
        evo = evolve(evo, tstep)
        push!(psis, evo.psi)
    end
    return psis
end
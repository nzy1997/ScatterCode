function time_evo(psi::ConstParticalNumberFermionState{T,PN}, H::CSCSimpleFermionHamiltonian, ttotal::Float64,tstep::Float64) where {T,PN}
    psis = ConstParticalNumberFermionState{T,PN}[]
    push!(psis, psi)
    for t in 0:tstep:ttotal
        psi = apply_hamiltonian(H, psi)
        push!(psis, psi)
    end
    return expm(-im * H * t) * ψ
end
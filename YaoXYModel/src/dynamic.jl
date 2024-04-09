# The XY model with driving term on the first site
# - `n` number of sites
# - `Jxy` coupling strength
# - `h` magnetic field strength
# - `Ω` driving strength
# - `frequency` frequency of the driving term
function drived_xymodel(n::Int, Jxy::Real, h::Real, Ω::Real, frequency::Real)
    H0 = Jxy/4 * sum([kron(n, i=>X, i+1=>X) + kron(n, i=>Y, i+1=>Y) for i=1:n-1]) +
        h/2 * sum([put(n, i=>Z) for i=1:n])
    H1 = Ω * put(n, 1=>X)/2
    return DrivedHamiltonian(H0, H1, frequency)
end

# get the Hamiltonian at time t
get_hamiltonian(dh::DrivedHamiltonian, t::Real) = dh.H0 + dh.H1 * (cos(dh.frequency * t))

# time evolve the state by Hamiltonian
# - `reg`: the quantum state,
# - `dh`: a DrivedHamiltonian instance,
# - `dt`: time step
# - `nstep`: number of steps
# - `operators`: a list of operators to measure
function simulate!(reg::AbstractRegister, dh::DrivedHamiltonian, dt, nstep::Int, operators)
    results = Vector{Vector{ComplexF64}}(undef, nstep+1)
    results[1] = Yao.expect.(operators, Ref(reg))
    for i = 1:nstep
        b = time_evolve(get_hamiltonian(dh, (i-0.5) * dt), dt)
        apply!(reg, b)
        results[i+1] = Yao.expect.(operators, Ref(reg))
    end
    return reg, results
end

# Fermionic annihilator (Jordan-Wigner transformation)
function fermion_an(n::Int)
    return [kron(n,[j=>Z for j = 1:i-1]...,i=>(X-im*Y)) for i in 1:n]/2
end

# Fermionic momentum operator
# ```math
# \eta_k = \sum_{j=1}^n U_{jk}^* c_j
# ```
function momentum_an(n::Int,U::Matrix)
    c=fermion_an(n)
    return [sum([U[j,k]'*c[j] for j in 1:n]) for k in 1:n]
end

# simulate the model with Yao
# - `n` number of sites
# - `k` target mode
# - `Ω` driving strength
# - `time_stop` total time
function sim_model(;
        n,
        k,
        Ω,
        time_stop,
        Jxy = 1.0,
        h = 1.0,
        time_step=0.1,
    )
    @assert 0 < k <= n "target mode k should be in 1 to $n, got: $k"

    reg = product_state(Yao.BitBasis.bit_literal(fill(1, n)...))
    model = XYModel(Jxy, h, n)
    A = fermionic_quadratic(model)
    eigs,vecs = eigen(A)
    frequency = -eigs[k]
    dh = drived_xymodel(n, Jxy, h, Ω, frequency)
    η = momentum_an(n,vecs)

    reg,dis = simulate!(reg, dh, time_step, round(Int,time_stop/time_step), vcat([η[i]'*η[i] for i in 1:n], [put(n,j=>Z) for j in 1:n]))
    return real.(hcat(dis...))'
end
"""
	continued_fraction(ϕ, niter::Int) -> Rational

obtain `s` and `r` from `ϕ` that satisfies `|s/r - ϕ| ≦ 1/2r²`
"""
continued_fraction(ϕ, niter::Int) = niter == 0 || isinteger(ϕ) ? floor(Int, ϕ) : floor(Int, ϕ) + 1 // continued_fraction(1 / mod(ϕ, 1), niter - 1)
continued_fraction(ϕ::Rational, niter::Int) = niter == 0 || ϕ.den == 1 ? floor(Int, ϕ) : floor(Int, ϕ) + 1 // continued_fraction(1 / mod(ϕ, 1), niter - 1)

function fraction_approximate(ϕ::Number, maxiter::Int; atol = 1e-4)
	for i in 1:maxiter
		fr = continued_fraction(ϕ, i)
		if abs(fr - ϕ) < atol
			return fr
		end
	end
	return continued_fraction(ϕ, maxiter)
end

function unitary_decomposition(u::AbstractMatrix)
	θ = 2 * acos(abs(tr(u)) / 2)
	if cos(θ / 2) < 1e-10
		θ = pi
		u1 = -(u) / im
		nx = abs(tr(u1 * [0 1; 1 0])) / 2
		expα = tr(u1 * [0 1; 1 0]) / (2 * nx)
		u1 = u1 / expα
		ny = real(tr(u1 * [0 -im; im 0])) / 2
		nz = real(tr(u1 * [1 0; 0 -1])) / 2
		return expα, θ, nx, ny, nz
	else
		expα = tr(u) / (2 * cos(θ / 2))
	end
	if abs(sin(θ / 2)) < 1e-10
		return expα, θ, 0.0, 0.0, 0.0
	else
		u1 = -(u / expα - cos(θ / 2) * I) / sin(θ / 2) / im
		nx = real(tr(u1 * [0 1; 1 0])) / 2
		ny = real(tr(u1 * [0 -im; im 0])) / 2
		nz = real(tr(u1 * [1 0; 0 -1])) / 2
		return expα, θ, nx, ny, nz
	end
end

function universal_check(ulist::Vector{Matrix{ComplexF64}}) #to be modified
	nmat = [0.0, 0.0, 0.0]
	for u in ulist
		expα, θ, nx, ny, nz = unitary_decomposition(u)
		cf = fraction_approximate(θ / pi, 10;atol = 1e-15)
		if cf isa Rational
			if cf.den > 10000
				nmat = hcat(nmat, [nx, ny, nz])
			end
		end
	end
	return rank(nmat) >= 2
end

function enlarge_ulist(ulist::Vector{Matrix{ComplexF64}})
	ulist_copy = copy(ulist)
	for u in ulist
		for u2 in ulist
			push!(ulist_copy, u * u2)
		end
	end
	return ulist_copy
end

function enlarge_ulist(ulist::Vector{Matrix{ComplexF64}}, maxiter::Int)
	ulist_copy = copy(ulist)
	for i in 1:maxiter-1
		ulist_copy = enlarge_ulist(ulist_copy)
	end
	return ulist_copy
end

function universal_check(ulist::Vector{Matrix{ComplexF64}}, enlarge_maxiter::Int)
	ulist_copy = enlarge_ulist(ulist, enlarge_maxiter)
	return universal_check(ulist_copy)
end

function random_unitary(n::Int)
	u, r = qr(randn((n, n)) + im * randn((n, n)))
	return Matrix(u)
end

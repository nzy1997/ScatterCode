using Test
using YaoXYModel: XYModel,fermionic_quadratic,drived_xymodel,momentum_an
using YaoXYModel.Yao
using LinearAlgebra

@testset "momentum_an" begin
    n=2
    Ω = 0.5
    k = 1
    Jxy = 1.0
    h = 1.0
    model = XYModel(Jxy, h, n)
    A = fermionic_quadratic(model)
    H = drived_xymodel(n, Jxy, h, Ω, -1.0)
    eigs,vecs = eigen(A)
    η = momentum_an(n,vecs)
    reg0 = product_state(Yao.BitBasis.bit_literal(fill(1, n)...))
    reg = Yao.apply(reg0, η[k]')
    @test Yao.expect(η[k]'*η[k],reg) ≈ 1
    E0 = Yao.expect(H.H0, reg0)
    E0_expected = eigen(Matrix(H.H0)).values[1]
    @test E0 ≈ E0_expected
    @test Yao.expect(H.H0, reg) ≈ eigs[k] + E0
end

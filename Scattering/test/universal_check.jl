using Test
using XYModelwithYao
using LinearAlgebra

@testset "fraction_approximate" begin
    @test fraction_approximate(1/3, 5) ≈ 1/3
    @test isapprox(fraction_approximate(12541/4333, 10) , 12541/4333; atol=1e-4)
    @show fraction_approximate(1/3+1e-10, 10)
end

@testset "unitary_decomposition" begin
    u=random_unitary(2)
    expα, θ,nx,ny,nz = unitary_decomposition(u)
    uapp = expα*(cos(θ/2)*I-im*sin(θ/2)*[nz nx-im*ny;nx+im*ny -nz])
    @test u ≈ uapp
end

@testset "enlarge_ulist" begin
    ulist = [random_unitary(2) for i in 1:2]
    @test (XYModelwithYao.enlarge_ulist(ulist,3) |> length) == 42
end

@testset "universal_check" begin
    ulist = [random_unitary(2) for i in 1:2]
    @show universal_check(ulist,2)
end

@testset "universal_check hardmard T" begin
    h=ComplexF64.(1/sqrt(2)*[1 1;1 -1])
    x=ComplexF64.([0 1;1 0])
    y=ComplexF64.([0 -im;im 0])
    z=ComplexF64.([1 0;0 -1])
    T=ComplexF64.([1 0;0 exp(pi/4*im)])
    ulist = [h,T]
    @test universal_check(ulist,4)==true

    ulist = [h,x,y,z]
    @test universal_check(ulist,3)==false
end
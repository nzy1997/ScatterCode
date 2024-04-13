using Test
using FermionSimulation

@testset "all_bases" begin
    sc = SpaceConfig((["↑", "↓"], ["1", "2", "3", "4"]))
    bases = all_bases(sc,3)
    @show bases
end

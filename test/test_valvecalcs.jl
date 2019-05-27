include("../src/valvecalculations.jl")

@testset "Dome pressures" begin
    #testing with R-value of zero to ignore the adjustment to PTRO
    @test domepressure_downhole(600, 0, 180) ≈ 750 atol = 5
    @test domepressure_downhole(727, 0, 120) ≈ 820 atol = 5
end

#TODO: full test of table calcs on an example well 

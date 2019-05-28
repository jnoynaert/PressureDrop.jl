include("../src/valvecalculations.jl")

@testset "Dome pressures" begin
    #testing with R-value of zero to ignore the adjustment to PTRO
    @test domepressure_downhole(600-14.65, 0, 180) ≈ 750 atol = 5 #function takes psig but test values are psia
    @test domepressure_downhole(727-14.65, 0, 120) ≈ 820 atol = 5
end


@testset "Thornhill-Craver" begin
    @test ThornhillCraver_gaspassage_simplified(900+14.7, 1100+14.7, 140, 16) ≈ 1127 atol=1
    @test ThornhillCraver_gaspassage_simplified(1200, 1100, 140, 16) == 0

    seats = [7*4, 32, 5*8, 3*16]
    Hernandez_results = [3.15, 4.11, 6.36, 9.39] .* 1000
    results = [0.827 * ThornhillCraver_gaspassage(420, 850, 150, s, 0.7) for s in seats]
    @test all(abs.(Hernandez_results .- results) .<= Hernandez_results .* 0.02) #2% tolerance to account for changing C_ds that aren't easily available (test cases use Winkler C_ds)

    @test ThornhillCraver_gaspassage(1200, 1100, 140, 16, 0.7) == 0
end

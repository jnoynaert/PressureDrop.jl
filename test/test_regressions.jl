
testpath = joinpath(dirname(dirname(pathof(PressureDrop))), "test/testdata")

@testset "Mayes 3-16 errors" begin
    default_GOR = 500
    max_depth = 8500

    GLpath = joinpath(testpath, "Mayes_3_16/Test_GLdesign_Mayes_3_16.csv")
    surveypath = joinpath(testpath, "Mayes_3_16/Test_survey_Mayes_3_16.csv")

    well = read_survey(path = surveypath, skiplines = 4, maxdepth = max_depth) #include IDs so that tailpipe and liner are properly accounted
    valves = read_valves(path = GLpath, skiplines = 1)

    model = WellModel(wellbore = well, roughness = 0.001,
                    valves = valves, pressurecorrelation = HagedornAndBrown,
                    WHP = 0, CHP = 0, dp_est = 25, temperature_method = "Shiu",
                    BHT = 160, geothermal_gradient = 0.8,
                    q_o = 0, q_w = 500,
                    GLR = 0, naturalGLR = 0,
                    APIoil = 38.2, sg_water = 1.05, sg_gas = 0.65)

    ## 2018-02-14 (domain error in H&B flow pattern)
    #duplicating the exact leadup from the loop where problem originally occurred, to avoid masking any problematic calculation
    oil = 41.3966666
    gas =  581.0
    water = 136.8298755
    TP = 138.0
    CP = 473.0
    gasinj = 533.1899861698156

    model.WHP = TP
    model.CHP = CP
    model.q_o = oil
    model.q_w = water
    model.naturalGLR = model.q_o + model.q_w <= 0 ? 0 : max( gas * 1000 / (model.q_o + model.q_w) , default_GOR * model.q_o / (model.q_o + model.q_w) )
    model.GLR = model.q_o + model.q_w <= 0 ? 0 : max( (gas + gasinj) * 1000 / (model.q_o + model.q_w) , model.naturalGLR)

    @test (gaslift_model!(model, find_injectionpoint = true, dp_min = 100) |> x-> x[1][end]) ≈ 678 atol = 1

    #2018-02-12 (no output)
    oil = 45.0
    gas =  584.0
    water = 149.0972222
    TP = 137.0
    CP = 474.0
    gasinj = 518.08798656154

    model.WHP = TP
    model.CHP = CP
    model.q_o = oil
    model.q_w = water
    model.naturalGLR = model.q_o + model.q_w <= 0 ? 0 : max( gas * 1000 / (model.q_o + model.q_w) , default_GOR * model.q_o / (model.q_o + model.q_w) )
    model.GLR = model.q_o + model.q_w <= 0 ? 0 : max( (gas + gasinj) * 1000 / (model.q_o + model.q_w) , model.naturalGLR)

    @test (gaslift_model!(model, find_injectionpoint = true, dp_min = 100) |> x-> x[1][end]) ≈ 697 atol = 1

    # nonconverging segment at 6395'
    oil = 2
    gas =  584.0
    water = 149.0972222
    TP = 137.0
    CP = 474.0
    gasinj = 518.08798656154

    model.WHP = TP
    model.CHP = CP
    model.q_o = oil
    model.q_w = water
    model.naturalGLR = model.q_o + model.q_w <= 0 ? 0 : max( gas * 1000 / (model.q_o + model.q_w) , default_GOR * model.q_o / (model.q_o + model.q_w) )
    model.GLR = model.q_o + model.q_w <= 0 ? 0 : max( (gas + gasinj) * 1000 / (model.q_o + model.q_w) , model.naturalGLR)

    @test (gaslift_model!(model, find_injectionpoint = true, dp_min = 100) |> x-> x[1][end]) ≈ 674 atol = 1
end
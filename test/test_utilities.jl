include("../src/utilities.jl")

@testset begin "Interpolation"
    md = [1.,2.,3.,4.,5.]
    inc = [1.,2.,3.,4.,1.]
    tvd = [1.,2.,3.,4.,5.]
    id = [1.,1.,1.,2.,1.]

    well = Wellbore(md, inc, tvd, id)

    property1 = [1.,1.,2.,3.,4.,1.]
    property2 = [1.,1.,1.,1.,2.,1.]

    @test interpolate(well, property1, 4.5) == 2.5

    points = [1.5, 4.5]
    @test interpolate_all(well, [property1, property2], points) == hcat([1.5, 2.5], [1, 1.5])
end

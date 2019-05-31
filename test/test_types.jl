@testset "Wellbore object" begin

md_bad = [-1.,2.,3.,4.]
md_good = [1.,2.,3.,4.]
inc = [1.,2.,3.,4.]
tvd_bad = [-1.,2.,3.,4.]
tvd_good = [1.,2.,3.,4.]
id = [1.,1.,1.,1.]

#implicit test for adding the leading 0,0 survey point:
w = Wellbore(md_good, inc, tvd_good, id)

@test w.id[1] == w.id[2]

#implicit test for allowing negatives:
Wellbore(md_bad, inc, tvd_bad, id, true)

try
    Wellbore(md_bad, inc, tvd_good, id)
catch e
    @test e isa Exception
end

try
    Wellbore(md_good, inc, tvd_bad, id)
catch e
    @test e isa Exception
end

end #testset for Wellbore object

@testset "Wellbore with valves" begin

    md = [1.,2,3,4]
    tvd = [2.,4,5,6]
    inc = [0.,2,3,4]
    id = [1.,1,1,1]
    starting_length = length(md)

    valve_md = [1.5,3.5]
    PTRO = [1000, 900]
    R = [0.07,0.07]
    ports = [16,16]
    valves = GasliftValves(valve_md, PTRO, R, ports)

    well = Wellbore(md, inc, tvd, id, valves)

    @test all(length.([md,tvd,inc,id]) .== starting_length)
    @test well.md == [0,1,1.5,2,3,3.5,4]
    @test well.tvd == [0,2.,3,4,5,5.5,6]
    @test well.inc == [0,0.,1,2,3,3.5,4]
    @test well.id == [1,1.,1,1,1,1,1]

end

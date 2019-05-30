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

    @test false #TODO: test adding valves to a new wellbore object

end

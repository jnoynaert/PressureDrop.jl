# A (sloppy) example of some actual analysis done using PressureDrop.jl to generate normalized pressure plots for several wells through time.


wells = ["<wellname>",  "<wellname>"]
baseline_pressure = Dict("<wellname>" => 3125, "<wellname>" => 3191)
patternplot = true
max_depth = 8500
start_date = "21-Jul-2019"

import Printf
using DataFrames
using RCall

function initialize_R()

    R"""
    library(DBI)
    library(dplyr)
    library(tidyr)

    print(paste0('Connecting to ', $db))

    connection <- DBI::dbConnect(odbc::odbc(), driver = 'SQL Server',
                                server = '<Aries Server>', database = '<Aries DB>')
    """
end


# note that variable interpolation in @R_str is actual variable interpolation, not string interpolation
function get_production(well)

    R"""

    query <- paste0("
    select d.D_DATE as Date, p.PropNum,
    isnull(case when d.Oil < 0 then 0 else d.Oil end,0) as Oil,
    isnull(case when d.Gas < 0 then 0 else d.Gas end,0) as Gas,
    isnull(case when d.Water < 0 then 0 else d.Water end,0) as Water,
    isnull(d.Press_FTP,0) as TP,
    isnull(d.Press_Csg,0) as CP,
    isnull(d.Press_Inj,0) as GasInj,
    d.REMARKS as Comments

    from ARIES.AC_DAILY as d
    inner join ARIES.AC_PROPERTY as p
    on d.PROPNUM = p.PROPNUM

    WHERE p.AREA = 'oklahoma'
    and p.LEASE like '%",$well,"%'
    and d.D_DATE >= '",$start_date,"'
    ORDER BY D_DATE asc")

    prod <- dbGetQuery(connection, query) %>%
                tidyr::fill(Water, TP, CP, GasInj, .direction = 'down') %>%
                replace_na(list(Oil = 0, Gas = 0))

    """

    @rget prod
    @assert length(unique(prod.PropNum)) == 1 "Query for $well accidentally pulled in $(length(unique(prod.PropNum))) propnums"

    return prod
end


function finish_R()

    R"""
    dbDisconnect(connection)
    """
end


using PressureDrop
using Gadfly


function calculate_pressures(well, prod, default_GOR = 500)

    well_path = read_survey(path = "Surveys/$well.csv", skiplines = 4, maxdepth = max_depth)
    valves = read_valves(path = "GL designs/$well.csv", skiplines = 1)

    model = WellModel(wellbore = well_path, roughness = 0.001,
                    valves = valves, pressurecorrelation = HagedornAndBrown,
                    WHP = 0, CHP = 0, dp_est = 25, temperature_method = "Shiu",
                    BHT = 160, geothermal_gradient = 0.8,
                    q_o = 0, q_w = 500,
                    GLR = 0, naturalGLR = 0,
                    APIoil = 38.2, sg_water = 1.05, sg_gas = 0.65)

    FBHP = Array{Float64, 1}(undef, nrow(prod))

    for day in 1:nrow(prod)
        if prod.Oil[day] + prod.Water[day] == 0
            FBHP[day] = NaN
        else
            model.WHP = prod.TP[day]
            model.CHP = prod.CP[day]
            model.q_o = prod.Oil[day]
            model.q_w = prod.Water[day]
            model.naturalGLR = model.q_o + model.q_w <= 0 ? 0 : max( prod.Gas[day] * 1000 / (model.q_o + model.q_w) , default_GOR * model.q_o / (model.q_o + model.q_w) ) #force minimum GLR
            model.GLR = model.q_o + model.q_w <= 0 ? 0 : max( (prod.Gas[day] + prod.GasInj[day]) * 1000 / (model.q_o + model.q_w) , model.naturalGLR)

            FBHP[day] = gaslift_model!(model, find_injectionpoint = true, dp_min = 100) |> x -> x[1][end]
        end
    end

    return model, FBHP, gaslift_model!(model, find_injectionpoint = true, dp_min = 100)
end

ticks = log10.([0.1:0.1:0.9;1.0:1:10] |> collect)
function plot_normP_data(production, FBHP, well)

    production = production[.!isnan.(FBHP), :]
    FBHP = FBHP[.!isnan.(FBHP)]

    initial_pressure = haskey(baseline_pressure, well) ? baseline_pressure[well] : maximum(FBHP)
    println("Baseline pressure: $initial_pressure")

    norm_rate = (production.Oil .+ production.Water) ./ (initial_pressure .- FBHP)
    println("Initial PI: $(norm_rate[1])")
    println("Max PI: $(maximum(norm_rate))")
    println("Final PI: $(norm_rate[end])")

    return plot(layer(x = 1:length(norm_rate), y = norm_rate, Geom.path), 
                Scale.y_log10(minvalue = 10^-1., maxvalue = 10^1, labels = y-> Printf.@sprintf("%0.1f", 10^y)), Scale.x_log10,
                Guide.ylabel("Normalized Total Fluid (bpd / ΔP)"),
                Guide.xlabel("Normalized Time (days)"),
                Guide.yticks(ticks = ticks),
                Guide.title(well))
end


Base.CoreLogging.disable_logging(Base.CoreLogging.Warn) #prevent dumping all of the @infos from normal PressureDrop functions

initialize_R()
GL_plots = Dict()
pressure_plots = Dict()
for well in wells
    println("Evaluating $well...")

    production = get_production(well)[1:end-1, :] #last day is usually 0 prod
    model, FBHPs, GL_data = calculate_pressures(well, production); #GL data is tuple of TP, CP, valve data

    #generate plots 1 by 1
    GL_plots[well] = plot_gaslift(model.wellbore, GL_data[1], GL_data[2], model.temperatureprofile, GL_data[3], well)
    pressure_plots[well] = plot_normP_data(production, FBHPs, well)

    println("$well complete")
end
finish_R()

#set_default_plot_size(6inch,8inch)

using Compose

# pressure plots
if patternplot
    # lay out the wells in question according to gunbarrel view
    set_default_plot_size(12inch,5inch)
    gridstack(Union{Plot,Compose.Context}[pressure_plots["Mayes 1706 2-16MH"]  pressure_plots["Mayes 1706 4-16MH"] ]) |>
                                          SVG("pattern.svg")


else
    set_default_plot_size(6inch,8inch)

    for well in wells
        GL_plots[well] |> SVG("GL plots/$well GL.svg")
    end
end

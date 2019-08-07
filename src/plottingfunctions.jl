using .Gadfly
using Compose: compose, context


"""
`plot_pressure(well::Wellbore, pressures, ctitle = nothing)`

Plot pressure profile for a given wellbore using the pressure outputs from one of the pressure traverse functions.

See `traverse_topdown` and `pressure_and_temp`.
"""
function plot_pressure(well::Wellbore, pressures, ctitle = nothing)

        plot(x = pressures, y = well.md, Geom.path, Theme(default_color = "deepskyblue"),
                Scale.x_continuous(format = :plain),
                Guide.xlabel("Pressure (psia)"),
                Scale.y_continuous(format = :plain),
                Guide.ylabel("Measured Depth (ft)"),
                Guide.title(ctitle),
                Coord.cartesian(yflip = true))
end


"""
`plot_pressure(m::WellModel, pressures, ctitle = nothing)`

Plot pressure profile for a given wellbore using the pressure outputs from one of the pressure traverse functions.

The `wellbore` field must be defined in the passed WellModel.

See `traverse_topdown` and `pressure_and_temp`.
"""
function plot_pressure(m::WellModel, pressures, ctitle = nothing)

        plot_pressure(m.wellbore, pressures, ctitle)
end


"""
`function plot_pressures(well::Wellbore, tubing_pressures, casing_pressures, ctitle = nothing, valvedepths = [])`

Plot relevant gas lift pressures for a given wellbore and set of calculated pressures.

See `traverse_topdown`, `casing_traverse_topdown`, and `pressure_and_temp`.
"""
function plot_pressures(well::Wellbore, tubing_pressures, casing_pressures, ctitle = nothing, valvedepths = [])

        plot(layer(x = tubing_pressures, y = well.md, Geom.path, Theme(default_color = "deepskyblue")),
             layer(x = casing_pressures, y = well.md, Geom.path, Theme(default_color = "springgreen")),
             layer(yintercept = valvedepths, Geom.hline(color = "black", style = :dash)),
                Scale.x_continuous(format = :plain),
                Guide.xlabel("Pressure (psia)"),
                Scale.y_continuous(format = :plain),
                Guide.ylabel("Measured Depth (ft)"),
                Guide.title(ctitle),
                Coord.cartesian(yflip = true))
end


"""
`plot_pressures(m::WellModel, tubing_pressures, casing_pressures, ctitle = nothing)`

Plot relevant gas lift pressures for a given wellbore and set of calculated pressures.

The `wellbore` field must be defined in the passed WellModel, with the `valves` field optional.

See `traverse_topdown`, `casing_traverse_topdown`, and `pressure_and_temp`.
"""
function plot_pressures(m::WellModel, tubing_pressures, casing_pressures, ctitle = nothing)

        valvedepths = m.valves === missing ? [] : m.valves.md

        plot_pressures(m.wellbore, tubing_pressures, casing_pressures, ctitle, valvedepths)
end


"""
`plot_temperature(well::Wellbore, temps, ctitle = nothing)`

Plot temperature profile for a given wellbore using the pressure outputs from one of the pressure traverse functions.

See `linear_wellboretemp` and `Shiu_wellboretemp`.
"""
function plot_temperature(well::Wellbore, temps, ctitle = nothing)

        plot(x = temps, y = well.md, Geom.path, Theme(default_color = "red"),
                Scale.x_continuous(format = :plain),
                Guide.xlabel("Temperature (°F)"),
                Scale.y_continuous(format = :plain),
                Guide.ylabel("Measured Depth (ft)"),
                Guide.title(ctitle),
                Coord.cartesian(yflip = true))
end


"""
`plot_pressureandtemp(well::Wellbore, tubing_pressures, casing_pressures, temps, ctitle = nothing, valvedepths = [])`

Plot pressure & temperature profiles for a given wellbore using the pressure & temperature outputs from the pressure traverse & temperature functions.

See `traverse_topdown`,`pressure_and_temp`, `linear_wellboretemp`, `Shiu_wellboretemp`.
"""
function plot_pressureandtemp(well::Wellbore, tubing_pressures, casing_pressures, temps, ctitle = nothing, valvedepths = [])

        pressure = plot(layer(x = tubing_pressures, y = well.md, Geom.path, Theme(default_color = "deepskyblue")),
                        layer(x = casing_pressures, y = well.md, Geom.path, Theme(default_color = "mediumspringgreen")),
                        layer(yintercept = valvedepths, Geom.hline(color = "black", style = :dash)),
                Scale.x_continuous(format = :plain),
                Guide.xlabel("psia"),
                Scale.y_continuous(format = :plain),
                Guide.ylabel("Measured Depth (ft)"),
                Guide.title(ctitle),
                Coord.cartesian(yflip = true),
                Theme(plot_padding=[5mm, 0mm, 5mm, 5mm]))

        placeholdertitle = ctitle === nothing ? nothing : " "
        temp = plot(x = temps, y = well.md, Geom.path, Theme(default_color = "red"),
                Scale.x_continuous(format = :plain),
                Guide.xlabel("°F"),
                Scale.y_continuous(labels = nothing),
                Guide.yticks(label = false),
                Guide.ylabel(nothing),
                Guide.title(placeholdertitle),
                Coord.cartesian(yflip = true),
                Theme(default_color = "red", plot_padding=[5mm, 5mm, 5mm, 5mm]))

        hstack(compose(context(0, 0, 0.75, 1), render(pressure)),
                compose(context(0.75, 1, 0.25, 1), render(temp)))
end


"""
`plot_pressureandtemp(m::WellModel, tubing_pressures, casing_pressures, ctitle = nothing)`

Plot pressure & temperature profiles for a given wellbore using the pressure & temperature outputs from the pressure traverse & temperature functions.

The `wellbore` and `temperatureprofile` fields must be defined in the passed WellModel, with the `valves` field optional.

See `traverse_topdown`,`pressure_and_temp`, `linear_wellboretemp`, `Shiu_wellboretemp`.
"""
function plot_pressureandtemp(m::WellModel, tubing_pressures, casing_pressures, ctitle = nothing)

        valvedepths = m.valves === missing ? [] : m.valves.md

        plot_pressureandtemp(m.wellbore, tubing_pressures, casing_pressures, m.temperatureprofile, ctitle, valvedepths)
end


"""
`plot_gaslift(well::Wellbore, tubing_pressures, casing_pressures, temps, valvedata, ctitle = nothing)`

Plot pressure & temperature profiles along with valve depths and opening/closing pressures for a gas lift well.

Requires a valve table in the same format as returned by the `valve_calcs` function.

See `traverse_topdown`,`pressure_and_temp`, `linear_wellboretemp`, `Shiu_wellboretemp`, `valve_calcs`.
"""
function plot_gaslift(well::Wellbore, tubing_pressures, casing_pressures, temps, valvedata, ctitle = nothing)

        valvedepths = valvedata[:,2]

        pressure = plot(layer(x = [valvedata[:,12];valvedata[:,13]], y = [valvedepths;valvedepths], Geom.point, Theme(default_color = "mediumpurple3")), #PVC and PVO
                        layer(x = tubing_pressures, y = well.md, Geom.path, Theme(default_color = "deepskyblue")),
                        layer(x = casing_pressures, y = well.md, Geom.path, Theme(default_color = "mediumspringgreen")),
                        layer(yintercept = valvedepths, Geom.hline(color = "black", style = :dash)),

                Scale.x_continuous(format = :plain),
                Guide.xlabel("psia"),
                Scale.y_continuous(format = :plain),
                Guide.ylabel("Measured Depth (ft)"),
                Guide.title(ctitle),
                Coord.cartesian(yflip = true),
                Guide.manual_color_key(nothing,
                                       ["TP", "CP", "Valves", "PVO/PVC"],
                                       ["deepskyblue", "mediumspringgreen", "black", "mediumpurple3"]),
                Theme(plot_padding=[5mm, 0mm, 5mm, 5mm]))

        placeholdertitle = ctitle === nothing ? nothing : " " #generate a blank title to align the top of the plots if needed
        temp = plot(x = temps, y = well.md, Geom.path, Theme(default_color = "red"),
                Scale.x_continuous(format = :plain),
                Guide.xlabel("°F"),
                Scale.y_continuous(labels = nothing),
                Guide.yticks(label = false),
                Guide.ylabel(nothing),
                Guide.title(placeholdertitle),
                Coord.cartesian(yflip = true),
                Theme(default_color = "red", plot_padding=[5mm, 5mm, 5mm, 5mm]))

        hstack(compose(context(0, 0, 0.75, 1), render(pressure)),
                compose(context(0.75, 1, 0.25, 1), render(temp)))
end

"""
`plot_gaslift(m::WellModel, tubing_pressures, casing_pressures, valvedata, ctitle = nothing)`

Plot pressure & temperature profiles along with valve depths and opening/closing pressures for a gas lift well.

Requires a valve table in the same format as returned by the `valve_calcs` function. The passed WellModel must also have the `wellbore`, `temperatureprofile`, and `valves` fields defined.

See `traverse_topdown`,`pressure_and_temp`, `linear_wellboretemp`, `Shiu_wellboretemp`, `valve_calcs`.
"""
function plot_gaslift(m::WellModel, tubing_pressures, casing_pressures, valvedata, ctitle = nothing)

        plot_gaslift(m.wellbore, tubing_pressures, casing_pressures, m.temperatureprofile, valvedata, ctitle)
end

using .Gadfly
using Compose: compose, context

"""
plot_pressure(well::Wellbore, pressures, ctitle = nothing)

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
plot_gaslift_pressures(well::Wellbore, tubing_pressures, casing_pressures, ctitle = nothing)

Plot relevant gas lift pressures for a given wellbore and set of calculated pressures.

See `traverse_topdown`, `casing_traverse_topdown`, and `pressure_and_temp`.
"""
function plot_gaslift_pressures(well::Wellbore, tubing_pressures, casing_pressures, ctitle = nothing)

        plot(layer(x = tubing_pressures, y = well.md, Geom.path, Theme(default_color = "deepskyblue")),
             layer(x = casing_pressures, y = well.md, Geom.path, Theme(default_color = "springgreen")),
                Scale.x_continuous(format = :plain),
                Guide.xlabel("Pressure (psia)"),
                Scale.y_continuous(format = :plain),
                Guide.ylabel("Measured Depth (ft)"),
                Guide.title(ctitle),
                Coord.cartesian(yflip = true))
end


"""
plot_temperature(well::Wellbore, temps, ctitle = nothing)

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
plot_pressureandtemp(well::Wellbore, pressures, temps, ctitle = nothing)

Plot pressure & temperature profiles for a given wellbore using the pressure & temperature outputs from the pressure traverse & temperature functions.

See `traverse_topdown`,`pressure_and_temp`, `linear_wellboretemp`, `Shiu_wellboretemp`.
"""
function plot_pressureandtemp(well::Wellbore, pressures, temps, ctitle = nothing)

        pressure = plot(x = pressures, y = well.md, Geom.path, Theme(default_color = "deepskyblue"),
                Scale.x_continuous(format = :plain),
                Guide.xlabel("psia"),
                Scale.y_continuous(format = :plain),
                Guide.ylabel("Measured Depth (ft)"),
                Guide.title(ctitle),
                Coord.cartesian(yflip = true),
                Theme(plot_padding=[5mm, 0mm, 5mm, 5mm]))

        temp = plot(x = temps, y = well.md, Geom.path, Theme(default_color = "red"),
                Scale.x_continuous(format = :plain),
                Guide.xlabel("°F"),
                Scale.y_continuous(labels = nothing),
                Guide.yticks(label = false),
                Guide.ylabel(nothing),
                Coord.cartesian(yflip = true),
                Theme(default_color = "red", plot_padding=[5mm, 5mm, 5mm, 5mm]))

        hstack(compose(context(0, 0, 0.75, 1), render(pressure)),
                compose(context(0.75, 1, 0.25, 1), render(temp)))
end



"""
plot_pressureandtemp(well::Wellbore, pressures, temps, ctitle = nothing)

Plot pressure & temperature profiles for a given wellbore using the pressure & temperature outputs from the pressure traverse & temperature functions.

See `traverse_topdown`,`pressure_and_temp`, `linear_wellboretemp`, `Shiu_wellboretemp`.
"""
function plot_gaslift(well::Wellbore, tubing_pressures, casing_pressures, temps, valvedepths = [], ctitle = nothing)

        pressure = plot(layer(x = tubing_pressures, y = well.md, Geom.path, Theme(default_color = "deepskyblue")),
                        layer(x = casing_pressures, y = well.md, Geom.path, Theme(default_color = "springgreen")),
                        layer(yintercept = valvedepths, Geom.hline(color = "black")),
                Scale.x_continuous(format = :plain),
                Guide.xlabel("psia"),
                Scale.y_continuous(format = :plain),
                Guide.ylabel("Measured Depth (ft)"),
                Guide.title(ctitle),
                Coord.cartesian(yflip = true),
                Theme(plot_padding=[5mm, 0mm, 5mm, 5mm]))

        temp = plot(x = temps, y = well.md, Geom.path, Theme(default_color = "red"),
                Scale.x_continuous(format = :plain),
                Guide.xlabel("°F"),
                Scale.y_continuous(labels = nothing),
                Guide.yticks(label = false),
                Guide.ylabel(nothing),
                Coord.cartesian(yflip = true),
                Theme(default_color = "red", plot_padding=[5mm, 5mm, 5mm, 5mm]))

        hstack(compose(context(0, 0, 0.75, 1), render(pressure)),
                compose(context(0.75, 1, 0.25, 1), render(temp)))
end

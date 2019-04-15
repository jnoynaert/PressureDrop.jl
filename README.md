# PressureDrop.jl [![Build Status](https://travis-ci.org/jnoynaert/PressureDrop.jl.svg?branch=master)](https://travis-ci.org/jnoynaert/PressureDrop.jl)
Julia package for computing multiphase pressure profiles for gas lift optimization of oil &amp; gas wells.

Currently calculates outlet-referenced models for producing wells using non-coupled temperature gradients.

Note that all calculations and inputs are currently in U.S. field units.

# Installation

TODO: register package & add instructions

# Example usage

TODO: add notebook exmaples

# Supported correlations

- Beggs and Brill 1973, with Payne correction factors. Best for inclined pipe.
- Hagedorn and Brown 1965, with Griffith and Wallis bubble flow correction. Best for high water cuts.

Neither correlation accounts for oil-water phase slip.

# Performance

The pressure drop calculations converge quickly enough in most cases that special performance considerations do not need to be taken into account during interactive use.

For bulk calculations, note that as always with Julia, the best performance will be achieved by wrapping any calculations in a function, e.g. a `main()` block, to enable proper type inference by the compiler.

Plotting functions are loaded separately by `import`ing or `using` the integrated `PressurePlots` module, to allow avoiding the startup overhead of the `Gadfly` plotting dependency.

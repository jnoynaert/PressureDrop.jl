# PressureDrop.jl
Package for computing pressure profiles for gas lift optimization of oil &amp; gas wells.

**Note that all calculations are currently in U.S. field units.**

# Installation

TODO: register package
TODO: publish minimalist DockerHub container with just Julia, this package, and dependencies

# Example usage

TODO. Add as a notebook

# Supported correlations

TODO

# Notes

Note about applicable range of pvt correlations available

# Performance

The pressure drop calculations converge quickly enough in most cases that special performance considerations do not need to be taken into account during interactive use.

For bulk calculations, note that as always with Julia code, the best performance will be achieved by wrapping any calculations in a function, e.g. a `main()` block, to enable proper type inference by the compiler.

The majority of startup and runtime delays are driven by the `Gadfly` dependency, which while excellent, incurs a significant precompilation penalty.

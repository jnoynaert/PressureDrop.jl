---
title: 'PressureDrop.jl: Pressure traverses and gas lift analysis for oil & gas wells'
tags:
  - Julia
  - petroleum engineering
  - gas lift
  - nodal analysis
authors:
 - name: Jared M. Noynaert
   orcid: 0000-0002-2986-0376
   affiliation: "1"
affiliations:
 - name: Alta Mesa Resources
   index: 1
date: 6 August 2019
bibliography: paper.bib
---

# Summary

The calculation of multiphase pressure profiles using empirical correlations is a common task in petroleum engineering research and practice, enabling diagnostics and transient analyses that otherwise depend on directly measured bottomhole pressure. Unfortunately, most options to perform these calculations depend on commerical software, and many software suites handle bulk calculations poorly or not at all. In addition, these commercial solutions typically have poor support for finding and modifying gas lift injection points for repeatedly modelling wells with that type of artificial lift. 

``PressureDrop.jl`` is a Julia package for computing multiphase pressure profiles for gas lifted oil and gas wells, developed as an open-source alternative to feature subsets of nodal analysis or RTA software such as Prosper, Pipesim, or IHS Harmony. It currently calculates outlet-referenced models for producing wells using non-coupled temperature gradients using industry-standard pressure correlations: Beggs and Brill [@Beggs:1973] with the Payne [@Payne:1979] correction, and Hagedorn and Brown [@Brown:1977] with the Griffith bubble flow [@Griffith:1961] correction, as well as the Ramey and Shiu temperature correlations [@Ramey:1962; @Shiu:1980]. Output plots are generated using `Gadfly.jl` [@Jones:2019].

In addition to being open-source, ``PressureDrop.jl`` has several advantages over closed-source applications for its intended use cases: (1) it allows programmatic and scriptable use with native code, without having closed binaries reference limited configuration files; (2) it supports dynamic recalculation of injection points and temperature profiles through time; (3) it enables duplication and modification of models and scenarios, including dynamic generation of parameter ranges for sensitivity analysis and quantification of uncertainty; (4) PVT or pressure correlation options can be extended by adding functions in Julia code (or C, Python, or R); (5) it allows developing wellbore models from delimited input files or database records.

# References
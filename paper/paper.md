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

The calculation of multiphase pressure profiles using empirical correlations is a common task in petroleum engineering research and practice. Unfortunately, most options to perform these calculations, especially in bulk, depend on commerical software. In addition, even commercial solutions have poor support for finding and modifying gas lift injection points for modelling purposes. 

``PressureDrop.jl`` is a Julia package for computing multiphase pressure profiles for gas lifted oil and gas wells, developed as an open-source alternative to feature subsets of commercial nodal analysis or RTA software such as Prosper, Pipesim, or IHS Harmony. It currently calculates outlet-referenced models for producing wells using non-coupled temperature gradients using the industry-standard Beggs and Brill`[@author:2001]` and Hagedorn and Brown`[@author:2001]` methods.

In addition to being open-source, ``PressureDrop.jl`` has several advantages over closed-source applications for its intended use cases: (1) it allows programmatic and scriptable use with native code, with no binaries consuming configuration files or awkward keyword specifications; (2) it supports dynamic recalculation of injection points and temperature profiles through time; (3) it enables duplication and modification of models and scenarios, including dynamic generation of parameter ranges; (4) PVT or pressure correlation options can be extended by adding functions in Julia code (or C, Python, or R); (5) developing wellbore models from delimited input files or database records.

# References
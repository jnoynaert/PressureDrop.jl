# Changes

## v1.0.2
- Added benchmarks to test suite
- Improved performance of core calculations by 15-20%
- Minor documentation updates

## v1.0.1
- Added JOSS paper draft
- Added assertions to guardrail pulling in negative production rates
- Added Windows test builds on Travis
- Modified `read_survey` to allow passing a valve object to automatically add the associated MDs

## v1.0

### Feature & interface
- Added gravity-based casing pressure calculations & gas lift valve tables
- Added gas lift plots
- Added support for bubble point pressure (rather than assuming oil is always under BPP)
- Added docs with examples
- Converted all user-facing meta-functions to use **psig** instead of absolute pressure, to reduce overhead when using field data
- Added struct-based arguments: all arguments reside in a WellModel struct and can be easily re-used/modified without 20-argument function calls

### Fixes
- Corrected issue with Griffith & Wallis bubble flow correction for Hagedorn & Brown pressure drops
- Corrected edge case issue in superficial velocity calculations for 0 gas production
- Improved test coverage

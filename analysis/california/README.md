# California Proposition 99 TRSC application

This directory implements the frozen California design. Configuration files
are fixed before TRSC post-treatment computation. The primary analysis uses
the 1970--1974 donor-only calibration block, 1975--1985 fitting block and
1986--1988 validation block. Structural q=1 and q=2 analyses remain separate.

Run the numbered scripts in order from the repository root. The scripts write
machine-readable outputs and never use California's post-treatment treated
outcomes for candidate, split, radius or certificate selection.

In the present workspace the R runtime is unavailable. Package-independent
data, static-SCM and geometry stages are executed by the accompanying Python
runner, while TRSC-dependent stages retain an explicit `missing_R_runtime`
status until run in an R-capable environment.

For the complete package-native run:

```sh
R CMD INSTALL .
Rscript analysis/california/run_all.R
```

The raw Stata file is intentionally excluded from version control. Set
`TRSC_CALIFORNIA_DTA` to a local canonical `smoking.dta`, place it under
`analysis/california/data/raw/`, or install the CRAN `tidysynth` package and let
`01_retrieve_validate.R` load its bundled `smoking` panel. Generated processed
data and RDS fit objects are also excluded; the compact audited tables,
certificates, solver summaries, and publication figures are retained.

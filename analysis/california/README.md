# California Proposition 99 TRSC application

This directory implements the frozen California design. Configuration files
are fixed before TRSC post-treatment computation. The primary analysis uses
the 1970--1974 donor-only calibration block, 1975--1985 fitting block and
1986--1988 validation block. Structural q=1 and q=2 analyses remain separate.

Run the numbered scripts in order from the repository root. The scripts write
machine-readable outputs and never use California's post-treatment treated
outcomes for candidate, split, radius or certificate selection.

The committed `output/` directory contains retained audited tables, solver
summaries and figures. The [completed audit manifest](output/run_manifest_completed_audit.json)
and [execution provenance](output/execution_provenance.csv) describe their
recorded R execution, manual pre-run patches and remaining provenance
limitations. The exact executed package working-tree hash was not retained.
An earlier Python-only status is preserved as historical evidence in the
audit records; it is not the current status of the retained results.

These outputs are plug-in identification results. The application does not
supply the additional statistical and transfer bounds needed for a nominal
confidence-region interpretation. This public release preserves the outputs
without rerunning or replacing them.

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

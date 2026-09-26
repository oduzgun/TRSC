# California data provenance

The audited run used the canonical `smoking.dta` panel (SHA-256
`3a9b6ead5803a86b0f22adc3fb71fac5ed2a8eba6620f745660a520acf3d68b5`).
The raw file is not committed. Supply it with the `TRSC_CALIFORNIA_DTA`
environment variable, place it at `analysis/california/data/raw/smoking.dta`,
or install the CRAN `tidysynth` package so the workflow can load its bundled
`smoking` data. The script validates the 39-state, 1970--2000 panel before use.

Generated `processed/smoking_long.csv` preserves the source rows and variables.
`processed/cigsale_wide.csv` reshapes only the `cigsale` outcome into one row
per year, with California first and the 38 donors in alphabetical order.
No values are imputed, scaled, winsorized, or manually edited. Ordinary ADH
predictor variables remain in the long file for provenance but are not turned
into extra TRSC time observations.

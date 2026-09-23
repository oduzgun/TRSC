# TRSC simulation and benchmark harness

This directory implements the **research design** proposed by the TRSC paper. It
is deliberately outside the stable R-package namespace.

It is not a results directory and contains no empirical or Monte Carlo findings.

## Purpose

The harness is intended to discriminate between methods and expose failure modes,
not to manufacture favorable examples. The implemented scenario labels are:

- `constant`: fixed predictive relationship;
- `affine`: linearly evolving relationship, in the null class of the q = 2
  penalty;
- `curved`: smooth nonlinear evolution;
- `abrupt`: discontinuous stress case intentionally outside the positive
  squared-difference smoothness theory;
- `exact_duplication`: a donor is copied exactly, making raw donor
  representations nonunique;
- `near_redundancy`: one donor is a noisy near-copy of another;
- `weak_excitation`: the fitting design is nearly collinear while the
  independent reference distribution remains full-dimensional;
- `reference_mismatch`: the calibration reference distribution differs from
  later fitting/validation donor-state support.

The design therefore mirrors the paper's requested stress dimensions: exact and
near redundancy, constant/affine/curved paths, abrupt changes, weak local
excitation, and reference-law mismatch.

## Primary losses

The principal simulation losses are representation-aware:

1. integrated path profile loss under the known reference second moment;
2. integrated held-out prediction loss under an independent evaluation law;
3. independent-reference cutoff prediction/profile loss;
4. held-out cutoff prediction loss, which differs from the reference-law loss
   in the support-mismatch design;
5. training reconstruction error, reported separately.

Raw donor-share mean squared error is returned only when the scenario explicitly
marks the generating representation as identified. It is `NA` for the exact
duplication design.

## Benchmarks

The current matched-information benchmarks are:

- static simplex synthetic control using the same fitting outcomes;
- period-specific simplex fitting (`kappa = 0`);
- q = 1 Euclidean raw-share difference regularization;
- q = 2 Euclidean raw-share difference regularization.

All four benchmarks receive both path metrics and an endpoint-based cutoff
evaluation. The raw-share penalties are intentionally non-representation-invariant
comparators. They are research benchmarks, not package recommendations.

## Running one design

From the repository root after installing TRSC:

```r
library(TRSC)
source("research/simulations/R/scenarios.R")
source("research/simulations/R/metrics.R")
source("research/simulations/R/benchmarks.R")
source("research/simulations/R/run-study.R")

run <- trsc_sim_run_one(
  scenario_name = "exact_duplication",
  seed = 1,
  q = c(1, 2),
  kappa = c(0.03, 0.3, 3)
)

run$selected
run$fit_metrics
run$benchmark_metrics
```

## Monte Carlo grid

```r
results <- trsc_sim_run_grid(
  scenarios = c("constant", "affine", "curved", "exact_duplication"),
  seeds = 1:100,
  output_dir = "research/simulations/output"
)
```

The output directory is ignored by Git. Seeds, scenario labels and tuning grids
are retained in the returned objects. Failed Monte Carlo runs are retained in
`simulation_failures.csv` and in the returned `failures` table; the grid does
not silently discard them or abort after the first numerical failure. Do not
commit simulation results until the research design and reporting protocol are
frozen.

## Interpretation rules

- Do not infer causal effects: every design here is pre-treatment.
- Do not call abrupt-path performance a guarantee for a bounded-variation class;
  that class is not established by the current theory.
- Do not interpret low empirical Gram rank as evidence that omitted population
  directions are irrelevant.
- Do not interpret raw-share instability under redundancy as predictive
  instability.
- Use oracle tuning only if it is explicitly labeled as an oracle diagnostic and
  applied comparably across methods.
- Unsuccessful fits and numerical issues must be retained in eventual Monte
  Carlo reporting rather than silently discarded.

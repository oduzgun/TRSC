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
  later fitting/validation donor-state support;
- `interior`: the true donor mixture stays strictly inside the simplex,
  stressing the full-dimensional tangent regime;
- `persistent_errors`: target innovations follow a persistent AR(1) process,
  stressing covariance inflation outside the independent-noise benchmark;
- `separated_face`: a paper-matched intrinsic hull
  `conv{-u, u, v ± theta h_j}` with a stable exposed target edge and a
  positive outward normal margin;
- `vanishing_margin`: the same exposed-edge geometry with a deliberately
  small normal margin, stressing the localization boundary;
- `face_change`: an exact-representation path switching between two exposed
  edges, explicitly outside any generic face-change theorem.

The design therefore mirrors the paper's requested stress dimensions: exact and
near redundancy, constant/affine/curved paths, abrupt changes, weak local
excitation, reference-law mismatch, interior targets, separated faces,
vanishing normal margins, face changes, and persistent errors.

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
- q = 2 Euclidean raw-share difference regularization;
- a local Gaussian-kernel simplex varying-coefficient smoother, which estimates
  a separate locally weighted simplex fit at each fitting date;
- a coordinate L1 total-variation simplex path, included as a change-point /
  discontinuity-oriented comparator.

All six benchmarks receive both path metrics and an endpoint-based cutoff
evaluation. The raw-share and total-variation coordinate penalties are
intentionally non-representation-invariant
comparators. They are research benchmarks, not package recommendations.

## Structural and design diagnostics

Each run now records fitting-design diagnostics separately from predictive loss:

- donor Gram rank and minimum eigenvalue on the intact fitting block;
- contrast Gram rank and minimum eigenvalue after differencing against one donor;
- learned calibration-geometry rank;
- population/reference Gram rank and minimum eigenvalue.

These quantities are descriptive stress diagnostics. In particular, low
empirical rank is not treated as evidence that omitted population directions are
irrelevant.

The harness also compares two **different validation targets** on the same
frozen TRSC candidate library:

1. the paper's signed local-linear cutoff-targeted selector;
2. ordinary historical validation using unweighted mean squared error over the
   retained validation dates.

Their scores are not interpreted as estimates of the same risk target.

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
are retained in the returned objects. Selection comparisons, design
diagnostics, and TRSC solver diagnostics retain the originating scenario and
seed. Failed Monte Carlo runs are retained in `simulation_failures.csv` and in
the returned `failures` table; the grid does not silently discard them or
abort after the first numerical failure. Do not
commit simulation results until the research design and reporting protocol are
frozen.

## Monte Carlo reporting protocol

After a grid has completed, source `R/reporting.R` and construct a reporting
object:

```r
source("research/simulations/R/reporting.R")
report <- trsc_sim_mc_report(results)
```

The report contains:

- mean loss estimates with Monte Carlo standard errors and normal-approximation
  intervals for the simulation mean;
- empirical 2.5% / 97.5% replication quantiles, reported separately from the
  Monte Carlo interval;
- per-scenario failure rates with an accounting check that successful and failed
  runs sum to all attempted seeds;
- cutoff-targeted versus ordinary-historical selection frequencies;
- fitting-design and calibration-rank summaries;
- TRSC solver iteration, primal/dual residual and simplex-feasibility
  summaries;
- a reproducibility manifest containing package version, R version, platform,
  scenario names, seeds, run counts and session information.

`trsc_sim_write_report_bundle()` writes these tables and the manifest to a
chosen output directory. Report bundles and numerical findings remain
uncommitted until the simulation design and reporting protocol are frozen.

Monte Carlo intervals quantify simulation uncertainty for the reported
experiment; they are not confidence intervals for the statistical estimand and
do not certify theorems or causal identification.

## Interpretation rules

- Do not infer causal effects: every design here is pre-treatment.
- Do not call abrupt-path or face-change performance a guarantee for a
  bounded-variation or generic changing-face class; those classes are not
  established by the current theory.
- Treat `separated_face` as theorem-matched only at the level of its stated
  structural primitives; `vanishing_margin` intentionally weakens a key
  localization condition.
- Do not interpret low empirical Gram rank as evidence that omitted population
  directions are irrelevant.
- Do not interpret raw-share instability under redundancy as predictive
  instability.
- Use oracle tuning only if it is explicitly labeled as an oracle diagnostic and
  applied comparably across methods.
- Unsuccessful fits and numerical issues must be retained in eventual Monte
  Carlo reporting rather than silently discarded.

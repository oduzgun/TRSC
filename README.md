# TRSC — Temporally Regularized Synthetic Controls

TRSC is research software for **temporally regularized synthetic-control predictive relationships**.

The central design rule is:

> **Predictive profiles are primary; raw donor coordinates are secondary representations.**

When donors are redundant, several simplex weight vectors can represent the same predictive relationship. TRSC therefore does not treat singular predictive geometry or nonunique raw donor weights as failures.

## Development status

This repository is under active methodological development and is currently private.

The pre-treatment statistical theory distinguishes prediction-space geometry from arbitrary Euclidean geometry in donor coordinates. For donor state (X), the prediction metric is the uncentered second moment

[
\Gamma_P = E_P[XX^\top],
]

so a change in donor representation is measured through

[
(w-v)^\top \Gamma_P (w-v).
]

The estimator is intended to support temporal difference orders (q=1) and (q=2) without requiring (Gamma_P) to be invertible.

The stable implementation now includes the **pre-treatment convex fitting core**. It does not yet include the paper's frozen candidate-selection/readout procedure and does not claim to deliver a post-treatment counterfactual, ATT, or treatment effect.

## Current stable foundation

- `trsc_data()`: validate and retain the declared raw panel and integrated raw blocks.
- `trsc_split()`: construct expanding chronological split assignments.
- `trsc_geometry()`: construct an empirical uncentered predictive geometry without inverting it.
- `trsc_profile_coordinates()`: map donor coordinates to empirical quotient coordinates.
- `trsc_profile_distance()`: evaluate predictive distance under the empirical geometry.
- `trsc_profile_equivalent()`: test empirical profile equivalence at a declared numerical tolerance.
- internal common difference operators for (q=1) and (q=2).
- `trsc_fit()`: simplex-constrained pre-treatment fitting under a supplied predictive geometry.
- `predictive_path()`: fitted predictive-relationship object carrying one numerical lift plus supplied-geometry coordinates; learned-metric null directions are not presented as population identification.
- `weights()`: one solver-selected donor-coordinate representation, explicitly secondary.
- standard `fitted()`, `residuals()` and `summary()` methods for fitted objects.

Empirical rank, null dimension and representation diagnostics are **diagnostics only**. They do not verify population coercivity, calibration complexity, continuation assumptions or causal identification.

## Pipeline

The intended development pipeline is

```text
data preparation
    ↓
chronological sample construction
    ↓
preprocessing
    ↓
predictive geometry
    ↓
representation-invariant TRSC fit
    ↓
candidate / tuning selection
    ↓
predictive path
    ↓
diagnostics
    ↓
experimental post-treatment identification layer
```

The fitting core is implemented. Candidate/readout selection, cutoff-targeted validation and any post-treatment identification layer remain separate future stages.

## Representation invariance

For two donor representations (w) and (v), zero predictive distance means they are equivalent at the geometry's fitted resolution. Adding or duplicating a predictively redundant donor may change a numerical weight vector while leaving the predictive profile unchanged.

The package exposes raw weights only as solver-selected coordinates. The optimizer does not add an identity ridge or other coordinate-based tie-break to manufacture uniqueness. If a numerical tie-break is introduced later for reproducibility, it will be documented as a **computational convention**, not as an identification result.

## Post-treatment boundary

Post-treatment work remains experimental. The theory currently separates:

1. continuation of the predictive-profile path;
2. ambiguity in evaluating a profile at an observed post-treatment donor state;
3. assumptions linking equivalent donor lifts across post-treatment dates.

No stable high-level function named `trsc_counterfactual()`, `trsc_treatment_effect()` or `trsc_att()` is provided.

Experimental work belongs under `research/` until these identification choices are frozen.

## Legacy code

Earlier TVSC scripts for version-1 blocks/scaling/metrics and the version-2 `tvsc.vmetrics()` architecture are preserved under `research/legacy/`. They are retained for traceability but are not loaded by the TRSC package because the current theory makes predictive geometry, rather than the earlier predictor-metric abstraction, the central object.

## Validation

GitHub Actions run:

- the preserved standalone contracts for data preparation and chronological splitting;
- package-level tests for predictive geometry, redundancy, donor permutation, singular fitting, the TRSC optimizer and common (q=1/q=2) temporal-difference infrastructure;
- `R CMD check` after regenerating roxygen documentation.

The repository contains no public release and no empirical application.

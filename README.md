# TRSC — Temporally Regularized Synthetic Controls

TRSC is research software for **temporally regularized synthetic-control predictive relationships**.

The central design rule is:

> **Predictive profiles are primary; raw donor coordinates are secondary representations.**

When donors are redundant, several simplex weight vectors can represent the same predictive relationship. TRSC therefore does not treat singular predictive geometry or nonunique raw donor weights as failures.

## Development status

This repository is under active methodological development and is currently private.

The pre-treatment statistical theory distinguishes prediction-space geometry from arbitrary Euclidean geometry in donor coordinates. For donor state (X), the prediction metric is the uncentered second moment

$
\Gamma_P = E_P[XX^\top].
$

so a change in donor representation is measured through

$
(w-v)^\top \Gamma_P (w-v).
$

The estimator is intended to support temporal difference orders (q=1) and (q=2) without requiring (Gamma_P) to be invertible.

The stable implementation now includes the **deployable pre-treatment procedure**: explicit calibration/fitting/validation role construction, donor-only calibration geometry, the convex fitting core, frozen terminal readouts, cutoff-targeted candidate selection, and optional dependence-transfer diagnostics. It does not claim to deliver a post-treatment counterfactual, ATT, or treatment effect.

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
- `trsc_diagnostics()`: descriptive fit/geometry/solver diagnostics with explicit non-certification labels.
- `trsc_readout()` / `trsc_readout_library()`: frozen convex terminal tapers, including the endpoint and dyadic bandwidth library.
- `trsc_select()`: signed local-linear cutoff-targeted validation over a frozen candidate library, with no final refit.
- `trsc_sample_plan()`: explicit earlier calibration segment, intact fitting block, validation segment, deterministic thinning, and gaps.
- `trsc_calibration_geometry()`: donor-only learned geometry from retained calibration dates.
- `trsc_candidate_library()`: finite prespecified `(q, kappa, h)` library, fitting each `(q, kappa)` path once before validation.
- `trsc_workflow()`: end-to-end pre-treatment calibration → fit → freeze → validate orchestration.
- `trsc_mixing_bound()`: evaluate the paper's four-term coupling bound under a supplied beta-mixing coefficient function.
- `trsc_geometric_mixing_design()`: compute the paper's conservative thinning/gap choices under an assumed exponential beta-mixing envelope.

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

The deployable pre-treatment workflow is implemented through frozen cutoff selection. Dependence-transfer arithmetic is available as an assumption-driven diagnostic; the package does not estimate or certify mixing. Any post-treatment identification layer remains a separate future stage.

## Deployable pre-treatment workflow

The theorem-aligned workflow keeps the data roles explicit:

```text
earlier raw calibration segment
    ↓ deterministic thinning
retained donor calibration states
    ↓
learned predictive geometry
    ↓ gap g0
intact contiguous fitting block
    ↓
finite (q, kappa) fit library
    ↓ frozen terminal h-readouts
    ↓ gap g1
retained validation dates
    ↓ signed local-linear cutoff score
selected frozen cutoff relationship
```

The candidate specification, calibration rule, bandwidth library, and tie order are fixed before validation is scored. There is no rolling refit and no final full-sample refit in this workflow.

For a dependent panel, the software can calculate

$
\Delta_{\mathrm{mix}}
=
(N_C-1)\beta(b_C)
+\beta(g_0)
+\beta(g_1)
+(v-1)\beta(b_V).
$

or the explicit exponential-envelope spacing rule supplied by the paper. These calculations are conditional on user-supplied mixing assumptions; they are not empirical certification of beta-mixing.

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
- workflow tests for role separation, calibration provenance, validation leakage, frozen candidate construction, deterministic tie handling and no-refit semantics;
- dependence-transfer tests for the exact four-term coupling expression and the explicit exponential-mixing spacing rule;
- `R CMD check` after regenerating roxygen documentation.

The repository contains no public release and no empirical application.

# Frozen TRSC paper to code crosswalk

Authority: Gate 5C `TRSC_Unified.tex`, SHA-256
`e3b6e199a27f580d6b28c226a5fb1b25c580d2d9f71be60ee869562afdd66d92`.
This checksum was verified against the authoritative workspace manuscript on
September 25, 2026. It supersedes the earlier crosswalk checksum as a
provenance correction only; no statistical result or estimator was changed.
Baseline repository: `a0b228f5abb6f28f6f7f0aff2ea0b5775cac87e9`.
This is an implementation contract, not a revision of the statistical theory.

## Baseline audit (historical, before this branch)

The 53-file baseline has a coherent Estimator A QP and chronological pipeline,
41 test blocks, three executable vignettes, standalone preparation tests and
a research-only pre-treatment simulation/reporting harness. All package R,
tests, workflow, metadata, examples, vignettes and research interfaces were
inspected. `research/legacy` is an explicitly excluded historical archive;
it must not be sourced by the package. No repository AGENTS instructions exist.

The post-treatment directory contains only a research-boundary note. There
are no evaluation LPs, continuation SOCPs, auxiliary constraints, confidence
certificates or effect-functional APIs. Generated help files are absent.
The current CI exercises only the earlier pipeline.

One implementation defect was corrected: baseline `geometry$B` dropped positive
eigenvalues below a numerical tolerance while the fit retains them in Gamma.
Use the raw forward factor X/sqrt(N), retaining every supplied direction;
rank tolerance remains diagnostic only. This is a software repair, not a
new metric, ridge, full-rank assumption, or theory change.

## Completed result map

| Paper object/result | Public interface | Internal mechanism | Verification |
|---|---|---|---|
| Estimator A; Proposition 3.1 | `trsc_geometry`, `trsc_fit`, profile helpers | existing OSQP QP; forward metric factor | existing fit/invariance tests; weak-direction regression |
| split–fit–freeze–validate | data, sample plan, calibration geometry, candidate library, selection, workflow | retained chronological pipeline | existing leakage, tie, split, mixing tests |
| terminal triangular readout | `trsc_readout`, `trsc_boundary_readouts`, library | convex level readout and deterministic same-taper preceding endpoint | readout and certificate tests |
| Theorem 8.1; evaluation fibers | `trsc_evaluation_set`, width, point-identification diagnostic | LP with profile equality and all simplex lifts | analytic duplicate/null examples |
| Theorem 9.1; dynamic set | `trsc_counterfactual_set`, `trsc_support` | anchored q1/q2 SOCP, independent future lifts | analytic and grid support tests |
| minimum curvature; sensitivity frontier | `trsc_min_radius`, `trsc_sensitivity` | root-energy epigraph SOCP; repeated support | q1 zero, q2 affine and Rmin example |
| auxiliary refinement, common latent witness | `trsc_auxiliary_refine` | affine/SOC measurement constraints before projection | fiber elimination, nesting, inconsistent measurement |
| Theorem 11.1 | anchor/geometry certificates, confidence set/support | empirical anchor SOCP; R + 2^(q-1)sqrt(H)epsilon | same-witness inclusion; singular geometry; certificate errors |
| Corollary 11.2 | certificate propagation diagnostics | joint W-weighted level/velocity set | horizon and joint-radius tests |
| Proposition 11.3 and modulus | `trsc_deployment_modulus` | two-lift SOCP | exact ambiguity and weak-direction example |
| Proposition M.1 | confidence diagnostics | outer radius accounting | population/outer/sandwich numerical fixtures |
| Lemma M.2 | internal interpolation diagnostic in tests | direct finite differences | vector-valued deterministic inequality |
| conditional-mean effect translation | `trsc_effect_interval` | certified endpoint translation | sign, probability and estimand tests |
| numerical input contract | all certificate constructors | explicit source/type/event/failure metadata; no oracle substitution | missing/inconsistent/conditional-certificate tests |

## Executable verification map

All paths below are relative to `tests/testthat/` unless stated otherwise.
The theory remains authoritative; tests verify implementations and fixtures,
not theorems in arbitrary designs.

| Result | Solver/helper | Exact test file |
|---|---|---|
| Proposition 3.1, Estimator A | existing `trsc_fit` OSQP assembly; full `geometry$B` | `test-fit.R`, `test-fit-invariance.R`, `test-conic-infrastructure.R` |
| Theorem 8.1 | `.trsc_solve_lp`, `.trsc_add_equality` | `test-identification.R`: evaluation fibers, singleton and duplicates |
| Theorem 9.1 | `.trsc_dynamic_program`, `.trsc_continuation_map`, `.trsc_solve_socp` | `test-identification.R`: independent grid enumeration, monotonicity and hull bounds |
| q1 persistence | first-difference map with R=0 | `test-identification.R`: profile fixed, evaluations free |
| q2 affine/minimum radius | second-difference map; `trsc_min_radius` epigraph | `test-identification.R`: affine path, exact Rmin=1, infeasibility below |
| Auxiliary refinement | `.trsc_add_auxiliary` on the future witness | `test-identification.R`: narrowing and incompatible measurements |
| Theorem 11.1 | `.trsc_confidence_program`, `.trsc_coverage` | `test-inference.R`: true-witness retention, missing/diagnostic/dimension/conditioning errors |
| Corollary 11.2 | `trsc_anchor_propagation`, exact 2-by-2 W calculations | `test-inference.R`: joint level/velocity and horizon amplification |
| Proposition 11.3 | `trsc_deployment_modulus`, two-lift conic program | `test-inference.R`: exact and near-null ambiguity; lower bound is theoretical, not estimated by the function |
| Proposition M.1 | `trsc_diagnostics` outer-radius accounting | `test-stress.R`, `test-post-diagnostics.R`: sandwich fixtures and allowances |
| Lemma M.2 | direct vector finite differences (no public estimator) | `test-stress.R`: stated inequality for k=2,...,20 |
| Common path-event conversion | `trsc_anchor_from_risk` | `test-stress.R`: Markov normalization and one shared readout event |
| Numerical upper endpoints | `.trsc_certify_upper`, `.trsc_primal_diagnostics` | `test-conic-infrastructure.R`: imperfect duals and failed-status fallback |
| Conditional-mean translation | `trsc_effect_interval` | `test-inference.R`: endpoint signs, additional failure allocation, rejected realized estimand |
| Simulation separation | observed-only runner versus truth-only metrics | `research/simulations/tests/test-post-smoke.R`: truth mutation and failure retention |

No frozen paper file is modified. Analytic finite-sample constants are supplied
through certificates, not estimated by the package. The crosswalk does not
claim automatic verification of assumptions, minimax bounds or population
metric recovery.

## Scientific boundaries

The fitting algorithm does not establish a causal bridge, a fitting-gap bound,
or a continuation law. Predictively selected q is not the structural q.
Population identified sets, plug-in sets and confidence regions have distinct
classes and printed descriptions. A finite collection of support evaluations
does not enumerate the whole convex set. No inferred ATT or realized missing
innovation is returned. Exact population geometry supplied by an oracle is
labelled as such in research work.

## Numerical contract

ECOS supplies approximate conic primal/dual candidates. A maximizing primal
value is descriptive only. Certified upper endpoints use weak duality with
an exact rational residual correction over simplex blocks (GMP arithmetic),
or the always-valid donor-hull bound if no usable dual is available. The
certificate uses the intended affine coefficients, not rounded solver
stationarity. Cone membership and final outward conversion are checked.
This is numerical certification conditional on the supplied mathematical
program and statistical certificates, not verification of causal assumptions.

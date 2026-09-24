# Frozen TRSC paper to code crosswalk

Authority: Gate 5C `TRSC_Unified.tex`, SHA-256
`8e73fc63116faede9a98b5e508af2ef05a630dda60c25abfd31773cc94743399`.
Baseline repository: `a0b228f5abb6f28f6f7f0aff2ea0b5775cac87e9`.
This is an implementation contract, not a revision of the statistical theory.

## Initial whole-repository audit

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

One implementation defect needs correction: `geometry$B` drops positive
eigenvalues below a numerical tolerance while the fit retains them in Gamma.
Use the raw forward factor X/sqrt(N), retaining every supplied direction;
rank tolerance remains diagnostic only. This is a software repair, not a
new metric, ridge, full-rank assumption, or theory change.

## Result map and implementation targets

| Paper object/result | Public interface | Internal mechanism | Verification |
|---|---|---|---|
| Estimator A; Proposition 3.1 | `trsc_geometry`, `trsc_fit`, profile helpers | existing OSQP QP; forward metric factor | existing fit/invariance tests; weak-direction regression |
| split–fit–freeze–validate | data, sample plan, calibration geometry, candidate library, selection, workflow | retained chronological pipeline | existing leakage, tie, split, mixing tests |
| terminal triangular readout | `trsc_readout`, library | convex readout; optional preceding endpoint | readout and certificate tests |
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

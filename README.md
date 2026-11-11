# TRSC — Temporally Regularized Synthetic Controls

TRSC estimates evolving **predictive relationships**. Predictive profiles are
primary; donor weights are coordinate representations that can be nonunique.

The package implements the frozen split–fit–freeze–validate procedure and
post-treatment functions for **assumption-indexed identified sets and honest
outer confidence regions**. These require explicit continuation and
causal-bridge assumptions; confidence regions additionally require numerical
error bounds and probability conditions. The package does not automatically
identify a causal effect.

## Public release and paper code

[Version 0.1.0](https://github.com/oduzgun/TRSC/releases/tag/v0.1.0) is the
first public software release. Its computational files are unchanged from
[`06e7c015f3579c69cdfb21569ea5689982524655`](https://github.com/oduzgun/TRSC/tree/06e7c015f3579c69cdfb21569ea5689982524655),
the code snapshot referenced by the manuscript. The release updates the
public documentation and citation metadata. Use the version tag or the
fixed commit when recording the code used in an analysis; `main` can change.

The repository contains the R package, examples, frozen simulation design
and runners, and the California application with retained audit outputs.
It is a **software release, not the complete paper replication archive**:
the full 73,400-attempt Monte Carlo outputs, separate 8,000-attempt inference
validation experiment, and post-freeze numerical recovery audit are not
included in this release. The California output manifest records the
provenance and limitations of its retained results.

## Installation

Install the released version from R:

```r
install.packages("remotes")  # once, if needed
remotes::install_github("oduzgun/TRSC@v0.1.0", build_vignettes = FALSE)
library(TRSC)
```

The package depends on `ECOSolveR`, `gmp`, `Matrix` and `osqp`; installing
dependencies from source may require a compiler and system libraries.
The research scripts and California files belong to the repository checkout,
rather than the installed package API.

## Where to start

| Task | Entry point |
|---|---|
| Learn the pre-treatment workflow | [Introduction](vignettes/trsc-introduction.Rmd) and [executable example](inst/examples/pre-treatment-workflow.R) |
| Compute post-treatment identified sets | [Guide](vignettes/post-treatment-identification.Rmd) and [example](inst/examples/post-treatment-identification.R) |
| Understand the inputs required for inference | [Inference guide](vignettes/honest-inference.Rmd) and [example](inst/examples/honest-inference.R) |
| Inspect the frozen simulation design | [Simulation guide](research/simulations/README.md), [cell registry](research/simulations/FROZEN_CELL_REGISTRY.csv) and [runner](research/simulations/run-frozen-study.R) |
| Inspect the California application | [Application guide](analysis/california/README.md) and [audited output manifest](analysis/california/output/run_manifest_completed_audit.json) |
| Cite or reuse the software | [Citation metadata](CITATION.cff) and [MIT license](LICENSE.md) |

## Workflow and API

| Layer | Functions |
|---|---|
| Prepare and split | `trsc_data()`, `trsc_split()`, `trsc_sample_plan()` |
| Prediction geometry | `trsc_geometry()`, `trsc_calibration_geometry()`, `trsc_profile_coordinates()`, `trsc_profile_distance()`, `trsc_profile_equivalent()` |
| Fit and freeze | `trsc_fit()`, `predictive_path()`, `trsc_readout()`, `trsc_boundary_readouts()`, `trsc_readout_library()`, `trsc_candidate_library()`, `trsc_select()`, `trsc_workflow()` |
| Evaluation identification | `trsc_evaluation_set()`, `trsc_evaluation_width()`, `trsc_point_identified()` |
| Dynamic identification | `trsc_counterfactual_set()`, `trsc_support()`, `trsc_bounds()`, `trsc_min_radius()`, `trsc_sensitivity()` |
| Structural measurements | `trsc_auxiliary_refine()` |
| Inference certificates | `trsc_anchor_certificate()`, `trsc_anchor_from_risk()`, `trsc_geometry_certificate()`, `trsc_anchor_propagation()` |
| Outer inference | `trsc_confidence_set()`, `trsc_confidence_support()`, `trsc_effect_interval()` |
| Limits and diagnostics | `trsc_deployment_modulus()`, `trsc_diagnostics()`, `trsc_mixing_bound()`, `trsc_geometric_mixing_design()` |

## Pre-treatment estimation

Estimator A uses an uncentered prediction metric
\(\Gamma_P=E[XX^\top]\) in the temporal penalty. The fitting QP uses OSQP;
there is no raw-weight ridge or inverse of the prediction metric. Practical
orders are q=1 and q=2.

The canonical workflow reserves an earlier donor-only calibration segment,
an intact fitting block, and validation dates. Candidate orders, penalties,
terminal triangular readouts and tie order are frozen before validation.
The local-linear validation score targets the cutoff. There is no final
full-sample refit. `trsc_split()` remains an earlier expanding-window utility;
`trsc_sample_plan()` supplies the canonical role separation.

For a q=2 boundary, `trsc_boundary_readouts(fit, h)` constructs the current
and preceding predictive readouts from the same frozen fit, using the same
triangular taper shifted back exactly one fitted date. It does not supply a
fitting-to-causal-cutoff transfer bound or validate structural q=2.

The forward factor `geometry$B = X/sqrt(N)` retains all supplied directions.
Numerical rank tolerance is diagnostic only. Its row dimension is N, rather
than a truncated estimated rank; clients should use profile helpers rather
than assume a minimal coordinate dimension.

## Post-treatment identification

```r
library(TRSC)
g <- trsc_geometry(matrix(c(0, 0, 1), 1,
                         dimnames = list(NULL, c("A", "B", "C"))))
x <- c(A = 0, B = 2, C = 1)
trsc_evaluation_set(g, x, representative = c(A = 1, B = 0, C = 0))

z <- trsc_counterfactual_set(g, rbind(x, x),
                            anchor = c(A = 1, B = 0, C = 0),
                            q = 1, R = 0)
trsc_bounds(z)
trsc_support(z, c(1, -1))
```

This is a **plug-in identified set**, not a confidence set. Calling
`type="population"` asserts that the supplied inputs really are population
geometry and causal boundary profiles; it does not estimate them.

Structural q and root-energy radius R are explicit. They are never inherited
from selected predictive q or inferred from kappa. q=1 at R=0 fixes the
predictive profile; evaluations can still vary among compatible lifts.
q=2 uses the preceding profile and cutoff level. Its affine path may leave
the hull, requiring a positive `trsc_min_radius()`. Below that radius,
a certified infeasible support problem returns an explicit status and -Inf.

An auxiliary restriction must have a shared-composition interpretation,
identified response maps and treatment exclusion (or a supplied deterministic
error bound). It intersects the **same latent path** before projection.
Intersecting separately computed outcome intervals is not equivalent.

## Honest inference and numerical certificates

An outer confidence set requires an anchor certificate and an **all-simplex-
chord** geometry certificate. Fitting-local calibration and empirical rank
diagnostics do not supply the latter. Anchor components separate estimation,
fitting-to-cutoff transfer, predictive-to-causal transfer and numerical error.
q=2 retains a joint W-weighted level–velocity region.

Certificate constructors check dimensions, provenance and probability
bookkeeping. They do not establish their supplied bounds. Unknown risk
envelopes, covariance/coercivity constants, transfer bounds and mixing
coefficients must be maintained numerical bounds or independently certified.
`diagnostic_only` inputs cannot support nominal coverage.

A common finite-library event allows validation selection without an extra
data split or independence assertion. It does not cover a later refit.
Default coverage is for fixed deployment experiments or unconditional random
deployment. Conditional-on-deployment requests require conditional certificates.
R remains structural sensitivity, distinct from the sampling enlargement.

ECOS produces numerical candidates for LP/SOCP problems. **An approximate
primal maximum is never reported as an honest upper endpoint.** Upper
enclosures use exact GMP rational arithmetic to validate cone duals and
correct their stationarity residual over simplex blocks. If no useful dual
is available, a donor-hull bound remains valid but can be wide. A numerical
certificate is conditional on the supplied mathematical program; it does not
validate statistical or causal assumptions. See
[the certification contract](inst/developer/NUMERICAL_CERTIFICATION.md).

`trsc_effect_interval()` translates a counterfactual confidence set using
a supplied treated **conditional mean functional**, or its confidence
interval. Observed treated realizations are not silently treated as means.

## Guides and examples

Seven executable vignettes cover the pre-treatment workflow, predictive
geometry, nonunique weights, post-treatment identification, continuation
orders, honest inference, and deployment amplification. Run:

```sh
Rscript inst/examples/pre-treatment-workflow.R
Rscript inst/examples/post-treatment-identification.R
Rscript inst/examples/honest-inference.R
```

## Verification and research materials

The test suite checks geometry, frozen selection, analytic identification
fixtures, same-witness outer inclusion, numerical certificates and twenty
software stress cases. GitHub Actions run standalone contracts, package
tests, examples, all vignettes, research smoke checks and R CMD check.

The computational base passed the
[R validation run for commit `06e7c015f357`](https://github.com/oduzgun/TRSC/actions/runs/36255724140),
including package tests, vignette builds, executable examples and
`R CMD check`. These checks concern software execution and do not rerun the
paper's full Monte Carlo study.

`research/simulations/` retains the existing thirteen pre-treatment scenarios,
six benchmark families and failure/reporting machinery. A separate
post-treatment extension separates observable inputs from truth, and a
[predeclared simulation design](research/simulations/SIMULATION_DESIGN.md)
covers identification, inference and out-of-theory stress. The runners and
frozen registries are included; the full Monte Carlo output archive is not.

`analysis/california/` contains the Proposition 99 analysis scripts,
configuration, audited tables, solver summaries and figures. Its regions
are plug-in identified sets, not sampling confidence regions. Read the
application's output manifest before interpreting the provenance of the
retained results.

The [paper-code crosswalk](inst/developer/paper-code-crosswalk.md) records
the frozen results, helpers and tests. Earlier development scripts remain
excluded under `research/legacy/`; they are not loaded by the package.

## Citation and license

Onur Düzgün (2026). *TRSC: Temporally Regularized Synthetic Controls*,
version 0.1.0. [GitHub release](https://github.com/oduzgun/TRSC/releases/tag/v0.1.0).
Machine-readable metadata are provided in [CITATION.cff](CITATION.cff).

The software is distributed under the [MIT license](LICENSE.md).

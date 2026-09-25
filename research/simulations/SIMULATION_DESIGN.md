# Frozen simulation design — package-completion gate

Design version 1; declared before large Monte Carlo execution. Only deterministic
fixtures and smoke repetitions may run in this gate. No performance findings
are asserted. Changes after execution begins must receive a dated amendment.

## Questions and classes

1. **Estimator:** Does frozen predictive TRSC recover evolving profiles, and
   when do redundancy, weak excitation or a cutoff gap defeat this goal?
2. **Identification:** When does one shared temporal budget tighten datewise
   bounds? How do persistence, affine continuation and positive minimum
   curvature differ? Which evaluation fibers survive exact profile knowledge?
3. **Inference:** How much do boundary uncertainty, slope error, future-chord
   calibration and numerical upper certification add to population widths?
4. **Stress/out-of-theory:** What fails under reference mismatch, invalid
   causal/gap certificates, abrupt continuation or amplified deployment?

The existing thirteen pre-treatment scenarios, six benchmark families and
failure-reporting machinery remain intact. Existing separated-face scenarios
are predictive pseudo-target experiments, not automatically causal bridges.
Post-treatment causal experiments are generated separately with an explicit
conditional-mean bridge and independent truth record.

## Planned parameter grid

Use baseline T=120, J=10, predictive rank r=3, sigma=.2, H=5; vary one axis at
a time, plus prespecified interactions: (rank,J), (gap,velocity,H),
(R-Rmin,amplification), (dependence,library size). This is a selected grid, not
the prohibitively large Cartesian product.

| Dimension | Planned levels |
|---|---|
| T | 60, 120, 240 |
| J | 3, 10, 40, 120 (including J>T) |
| Predictive rank | 1, min(3,J), J; exact and near redundancy |
| Near-duplicate distance | 0, 1e-6, 1e-3, .1 |
| Noise sigma | 0, .1, .5, 1 |
| Horizon H | 1, 5, 10, 20 |
| Final fitting-to-cutoff gap | 3, 10, 25 dates; feasible splits only |
| q | structural 1 and 2, reported separately from predictive selection |
| Boundary velocity | 0, .05/T, .2/T; precise/weak estimation through noise and readout |
| Continuation paths | persistence, affine, curved, abrupt, hull exit |
| R | 0, Rmin, Rmin+.01D, Rmin+.1D, Rmin+D, 4 sqrt(H) D |
| Deployment multiplier | 1, 10, 100, 1000 in weak/null directions |
| Candidate count | 1, 4, 12, 24 frozen (q,kappa,h) combinations |
| Dependence | independent; stationary refresh process rho=.5,.9 |

D denotes the maintained profile-diameter bound. Unsupported combinations
(e.g., too few dates for a split, an affine historical simplex path that exits
its hull, an impossible requested rank) are design errors retained as failures,
never silently modified. R<Rmin cases are intentional structural infeasibility
checks and are reported separately from optimization failures.

## DGPs and information boundaries

The post-treatment generator uses bounded unit-norm donor response columns F,
stationary sign states with identity second moment, and an explicit simplex
relationship path. Pre-treatment observations are X_t w_t plus independent
noise; post-treatment conditional means are x_h'w_h. Deployment is observed
and may leave the reference support. The true Gamma, profile path, boundary
level/velocity, weights, mean path, minimum radius and population support
oracle live only under `truth`. Observables contain pre-treatment panels,
role indices, and deployment X. No realized treated outcome is substituted for
a treated conditional mean.

The post-treatment IDs are persistence, affine, hull_exit, exact_ambiguity,
near_redundancy, weak_excitation, full_rank, curved and abrupt. Curved/abrupt
continuations are in-class only when the declared R covers their true energy;
otherwise label them stress/out-of-theory. High amplification alone is an
allowed deployment regime, not a violation of the conditional-mean bridge.

## Oracle, semi-feasible and fully feasible layers

* **Oracle:** true geometry and causal anchors; population identified supports.
* **Semi-feasible:** fitted frozen centers with known population geometry and
  independently supplied anchor envelopes. The true geometry is openly used.
* **Fully feasible under maintained bounds:** empirical calibration geometry,
  frozen fitted centers and a predeclared bounded-response protocol. Universal
  hull envelopes can be intentionally wide. Sharper supplied risk/gap/chord
  certificates may be compared only with their provenance and probability
  allocation recorded. No realized error or oracle roughness is substituted
  for an unavailable certificate.

Default smoke inference uses known design diameter bounds; it is a software
check of the conservative construction, not evidence of statistical efficiency.
The callable feasible runner accepts only `observed` and `protocol`, never
`truth`. A truth-mutation regression test guards this boundary. Calibration
geometry from the empirical sample is not mislabeled population truth.

## Benchmarks and metrics

Retain static simplex, period-specific simplex, raw-weight quadratic q1/q2,
local varying simplex, and total-variation families. Raw-weight benchmarks
are explicitly representation-dependent research comparisons, not variants
of the canonical package estimator. Retain ordinary historical validation
versus cutoff-targeted validation as a selection comparison.

Post-treatment comparisons: selected fixed lift; fixed-profile fibers;
datewise donor hulls; q1/q2 dynamic sets; oracle population inputs; plug-in
inputs; honest outer regions. Optional exact auxiliary restrictions are a
separate stronger-information regime. Do not conflate their assumptions.

Primary outputs:

* Integrated predictive loss and cutoff profile loss; held-out predictive loss.
* Population identified width, statistical enlargement and total confidence
  width, each annotated when only a numerical enclosure is available.
* Support-endpoint coverage for declared directions e_h, their negatives,
  all-ones cumulative and equal-weight mean directions.
* Simultaneous true-mean-path coverage of coordinate bands; this is not
  coverage of the full convex identified set.
* Whole-set coverage only with an analytic set comparison or a verified
  sufficient input event. A sufficient-event frequency is a lower bound on
  whole-set coverage, not its exact frequency. Otherwise record NA, not a
  grid-based substitute.
* Rmin lower/upper certificates; joint support versus the sum of datewise
  dynamic-coordinate supports and versus unrestricted hull support.
* Width/unrestricted-hull-width; amplification modulus; joint versus marginal
  level/velocity propagation; library probability cost.
* Solver status, raw residuals, certified dual correction, certified gap when
  available, fallback frequency, wall-clock time and total failure accounting.

A certified upper enclosure alone cannot determine whether an outer endpoint
covers the *exact* population support if numerical enclosures overlap. Store
the comparison as certified-covered / certified-not-covered / unresolved.
Never turn unresolved numerical status into statistical success.

## Reproducibility and planned reporting

Planned Monte Carlo: 1000 repetitions per primary cell, 200 for expensive
stress cells, seeds `1000000 + 10000*scenario_index + repetition`. Scenario
indices are fixed by the declared ordered ID lists in generators; cell IDs
include every parameter. Smoke seeds are 11, 12 and 101 and are excluded from
scientific summaries. Do not choose cells or repetitions after inspecting
coverage. Record every attempted repetition including errors and warnings.

Record Git SHA, dirty status, package version, R/session info, ECOSolveR/OSQP/GMP
versions, RNG kind, protocol and certificate provenance, and full parameter
list. Save results only to ignored results/output directories. For coverage
report binomial Monte Carlo uncertainty and denominators; for loss and width
retain existing mean/MCSE reporting. Failure rates use all attempts.

Planned figures: predictive loss versus rank/redundancy; dynamic tightening
versus R; q2 uncertainty versus horizon; population and confidence sensitivity
curves; width versus deployment modulus; coverage versus gap/noise/dependence;
numerical correction and runtime. Planned tables: benchmark losses; coverage
and widths by information layer; failures/infeasibility/fallbacks; certificate
inputs. No large-run results belong in this design document or CI.

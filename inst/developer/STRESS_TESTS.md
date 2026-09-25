# Software stress-test register

All cases below have executable checks. Wide sets and structural infeasibility
are scientific/numerical outcomes, not automatically software defects.

| # | Case | Expected classification | Evidence under tests/testthat |
|---|---|---|---|
| 1 | Exact duplicate donors | correct result (augmented invariance) | test-identification.R, test-fit-invariance.R |
| 2 | Near duplicates | correct result; weak direction retained | test-conic-infrastructure.R, test-inference.R |
| 3 | Singular calibration geometry | meaningfully wide set | test-stress.R, test-geometry.R |
| 4 | J>T | correct result | test-stress.R |
| 5 | Simplex-boundary anchor | correct result | test-stress.R, test-identification.R |
| 6 | q1 R=0 | correct result with evaluation fibers | test-identification.R |
| 7 | q2 R=Rmin | correct result | test-identification.R, test-stress.R |
| 8 | Affine q2 path outside hull | structural infeasibility below positive Rmin | test-identification.R |
| 9 | Very large R | correct donor-hull saturation | test-stress.R |
| 10 | Exact evaluation ambiguity | meaningfully wide set | test-identification.R, test-inference.R |
| 11 | Deployment amplification | meaningfully wide set | test-inference.R |
| 12 | Nonunique empirical minimizer | correct fit; no manufactured uniqueness | test-stress.R; certificate contract requires uniform envelope |
| 13 | Inaccurate/failed solver status | certified upper enclosure, possibly wide | test-conic-infrastructure.R |
| 14 | Noise-free design | correct result | test-stress.R |
| 15 | Donor scale heterogeneity | correct forward-seminorm/hull result | test-stress.R |
| 16 | Donor permutation | correct invariant support | test-identification.R, test-fit-invariance.R |
| 17 | Equivalent raw input ordering | correct identical fitted predictions | test-stress.R |
| 18 | Missing certificate | clean error | test-inference.R |
| 19 | Inconsistent certificate dimensions | clean error | test-inference.R |
| 20 | Unsupported conditional coverage | clean error | test-inference.R |

Minimum-radius and support values are numerical enclosures. An inconclusive
solver status must not be interpreted as structural feasibility. Redundant
lifts, evaluation-distinct lifts and genuinely different predictive minimizers
remain different cases; the program does not equate their numerical endpoints.

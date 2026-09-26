# Scientific audit of the completed California TRSC run

## Audit scope and provenance

The audit uses the frozen configuration and the completed outputs supplied in
`california_trsc_results.zip` (SHA-256
`c86dc8fba654cc96571cb5f0656c382390439ed98220c705925ac228c5999923`).
The package base commit recorded before execution is
`64b23b38cbc339492547a752759d994db501a7cc`. The stale pre-run manifest reports
working-tree checksum
`a696470d757cd88d0c9460a3278eb6a5b82593540caffb8e70ece00e0f01e68e`, but this
is **not** the checksum of the package that generated the completed results.
The terminal history shows that the macOS source tree was subsequently patched
to decode the labelled Stata state variable and to add/export
`trsc_boundary_readouts()`. That post-patch tree was not archived, so its exact
checksum is unavailable. The manuscript checksum is
`e3b6e199a27f580d6b28c226a5fb1b25c580d2d9f71be60ee869562afdd66d92`,
and the dataset checksum is
`3a9b6ead5803a86b0f22adc3fb71fac5ed2a8eba6620f745660a520acf3d68b5`.

The successful execution used TRSC 0.1.0 under an R 4.5 arm64 framework on
macOS, with ECOSolveR 0.6.1 and gmp 0.7-5.1. The exact R patch release and the
installed `osqp` version were not recorded. These details are now stored in
`output/execution_provenance.csv`.

A provenance inconsistency must be recorded. The archived `run_manifest.json`,
`all_endpoint_status.csv`, `fairness_audit.csv`, `figure_status.csv`, and four
placeholder tables still describe the earlier Python-only state. The numbered
stage-status files and the completed R-generated outputs postdate those files
and show stages 01--08 and 10 completed, with stage 09 completed without an
honest set. The stale files were preserved; a corrected audit manifest and
reconciliation table were added rather than rewriting the raw record.

## 1. Pre-treatment selection

The validation scores are:

### q=1

| kappa | h=1 | h=2 | h=4 |
|---:|---:|---:|---:|
| 0.03 | 3.378447 | 5.598421 | 8.580046 |
| 0.10 | 3.474496 | 5.954791 | 9.044203 |
| 0.30 | 3.499612 | 6.230636 | 9.321894 |
| 1 | 3.939098 | 6.832548 | 10.207313 |
| 3 | 3.496079 | 6.029035 | 9.625359 |
| 10 | 2.501695 | 3.932157 | 6.835249 |
| 30 | **2.452774** | 3.858338 | 7.688276 |
| 300 | 3.005799 | 4.808109 | 9.581541 |

The original grid ended at 30, so q=1 triggered the prespecified upper-edge
rule. The presence of the 300 row verifies that the grid was expanded before
the selection object was frozen for the post-treatment scripts. The selected
candidate is `(q=1,kappa=30,h=1)`, and 30 is interior to the expanded grid.

### q=2

| kappa | h=1 | h=2 | h=4 |
|---:|---:|---:|---:|
| 0.03 | 4.815425 | 6.067506 | 8.549012 |
| 0.10 | **4.541497** | 6.289317 | 8.996769 |
| 0.30 | 4.576839 | 6.692959 | 10.166158 |
| 1 | 5.778143 | 7.492266 | 11.489792 |
| 3 | 5.857637 | 7.055845 | 11.835603 |
| 10 | 6.416454 | 7.629952 | 12.702200 |
| 30 | 6.841135 | 8.264283 | 13.280949 |

The selected candidate is `(q=2,kappa=0.1,h=1)`, which is not on a grid edge.
The code and outputs keep q=1 and q=2 in separate libraries and output files;
the lower q=1 score was not used to eliminate q=2 as a structural analysis.

The fitting RMSEs are \(6.47\times10^{-6}\) for q=1 and
\(1.70\times10^{-6}\) for q=2. These near-zero values reflect the flexible,
rank-deficient design and must not be read as post-treatment identification.

## 2. q=2 boundary state

The selected q=2 readout has `h=1`, level endpoint 1985, and preceding endpoint
1984. The predictive metric coordinates are:

| coordinate | beta(-1) | beta(0) | velocity |
|---|---:|---:|---:|
| reference_1 | 45.620663 | 45.178329 | -0.442334 |
| reference_2 | 47.022824 | 46.572957 | -0.449867 |
| reference_3 | 49.035971 | 48.554485 | -0.481485 |
| reference_4 | 49.660363 | 49.149736 | -0.510628 |
| reference_5 | 50.884828 | 50.363673 | -0.521155 |

The velocity norm is 1.078058. The level and previous representative weight
vectors each sum to one. The largest level weights are Utah 0.2643, New Mexico
0.1022, Idaho 0.0323, South Dakota 0.0310, and North Dakota 0.0309. Those weights
are computational lifts; the five-coordinate profile is the primary object.

The package helper constructed both readouts from the same selected q=2 fit and
the same taper shifted by one fitted date. With `h=1`, this is exactly the 1985
endpoint and the 1984 endpoint. It does not resolve the 1985-to-1988
fitting-to-cutoff gap.

The terminal history confirms that this helper was manually added to the local
package immediately before the successful run. Although the exact helper source
was not included in the results archive, the saved boundary values equal the
1985 and 1984 columns of the selected fit's metric-coordinate and representative-
weight files. Thus the h=1 boundary used in this run is internally verifiable;
the unarchived source patch remains a reproducibility limitation for other
bandwidths.

## 3. Minimum continuation radii

For q=1, the returned interval is exactly `[0,0]` with status
`analytic_persistence`, as required by persistence feasibility.

For q=2, the certified interval is `[0,1.478804e-08]`; the numerical estimate
is `-1.16e-11`, effectively zero, and the solver status is optimal. Therefore
the affine q=2 continuation is feasible in the supplied plug-in geometry. The
application does not exhibit a substantively positive q=2 minimum curvature.

## 4. Dynamic identified sets

The q-specific radii are `R_min + c*S_q,H`. The multiplier grid ultimately is

`0, .25, .5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024`.

For q=1, `S=4.073408e-05`; for q=2, `S=0.002804297`. Extension beyond 16
follows the frozen doubling rule. Both `saturated` flags are false.

### Minimum-radius coordinate bounds

| year | q1 lower | q1 upper | q1 width | q2 lower | q2 upper | q2 width |
|---:|---:|---:|---:|---:|---:|---:|
| 1989 | 84.552 | 117.074 | 32.522 | 79.510 | 108.752 | 29.242 |
| 1990 | 80.790 | 114.146 | 33.356 | 74.639 | 104.225 | 29.586 |
| 1991 | 76.199 | 114.002 | 37.803 | 70.507 | 103.007 | 32.500 |
| 1992 | 75.480 | 116.902 | 41.422 | 69.517 | 105.047 | 35.530 |
| 1993 | 75.842 | 113.324 | 37.482 | 69.898 | 167.600 | 97.702* |
| 1994 | 74.538 | 110.979 | 36.440 | 68.932 | 103.361 | 34.429 |
| 1995 | 71.122 | 115.853 | 44.731 | 65.071 | 105.476 | 40.405 |
| 1996 | 71.098 | 113.353 | 42.255 | 54.000 | 101.013 | 47.013* |
| 1997 | 71.904 | 113.996 | 42.092 | 63.823 | 99.368 | 35.545 |
| 1998 | 66.212 | 113.952 | 47.741 | 54.754 | 95.608 | 40.854 |
| 1999 | 66.151 | 110.799 | 44.647 | 55.794 | 90.452 | 34.658 |
| 2000 | 60.869 | 105.283 | 44.413 | 50.681 | 84.595 | 33.914 |

`*` The q=2 1993 upper and 1996 lower use simplex-hull fallback certificates,
not sharp conic dual certificates.

At q=1's minimum radius, median width is 41.76, maximum width 47.74, and median
hull fraction 0.337. At the largest radius, median width is 43.81 and the
largest hull fraction is 0.416. q=1 therefore remains materially tighter than
the donor hull throughout the examined grid.

At q=2's minimum radius, median width is 35.09. The reported maximum 97.70 is
fallback-driven. At the largest radius, median width is 79.68; the horizon-12
width is 106.73 and its hull fraction is 0.924. q=2 approaches the hull at long
horizons but does not formally saturate.

Exact support functions must be monotone under nested radii. q=1 has three raw
enclosure deviations, all tiny (largest upper decrease 0.00176). q=2 has 29
raw deviations, several large, aligned with poor primal diagnostics and
fallback certificates. This is a numerical-enclosure issue, not evidence that
the mathematical identified sets are nonnested. The paper should use the
audited monotone-closure figure and retain the raw curves in the supplement.

## 5. Evaluation fibers

Primary-geometry widths are:

| year | q1 width | q2 width | donor-hull width |
|---:|---:|---:|---:|
| 1989 | 32.522 | 30.028 | 115.900 |
| 1990 | 33.342 | 30.936 | 129.100 |
| 1991 | 37.803 | 35.235 | 116.900 |
| 1992 | 41.422 | 38.796 | 112.600 |
| 1993 | 37.482 | 34.503 | 111.400 |
| 1994 | 36.440 | 32.566 | 114.300 |
| 1995 | 44.731 | 39.296 | 123.300 |
| 1996 | 42.153 | 37.008 | 125.000 |
| 1997 | 42.079 | 37.510 | 129.800 |
| 1998 | 47.739 | 42.464 | 131.500 |
| 1999 | 44.638 | 39.205 | 127.800 |
| 2000 | 44.409 | 39.296 | 115.500 |

No primary or full-pre diagnostic fiber is point identified at the package's
tolerance. For q=1 the primary min/median/max widths are
32.52/41.75/47.74; for q=2 they are 30.03/37.26/42.46. The full-pre diagnostic
min/median/max widths are 1.13/2.66/3.51 for q=1 and 2.50/6.83/9.24 for q=2.

Three distinct facts must not be conflated:

1. The calibration geometry has many nonunique donor representations.
2. Positive R permits nonunique future predictive-profile paths.
3. Positive fiber widths prove nonunique deployment evaluations at a fixed
   supplied profile.

The third conclusion comes from the endpoint programs, not from empirical rank
alone.

## 6. Deployment modulus

For the primary geometry, the per-year exact moduli at epsilon zero are:

| year | omega(0) |
|---:|---:|
| 1989 | 41.103 |
| 1990 | 43.592 |
| 1991 | 44.381 |
| 1992 | 51.742 |
| 1993 | 42.855 |
| 1994 | 46.166 |
| 1995 | 51.145 |
| 1996 | 50.102 |
| 1997 | 59.515 |
| 1998 | 54.806 |
| 1999 | 53.308 |
| 2000 | 56.378 |

Summary over deployment years:

| epsilon | primary median | primary max | full-pre diagnostic median | full-pre diagnostic max |
|---:|---:|---:|---:|---:|
| 0 | 50.624 | 59.515 | 17.839 | 20.925 |
| 0.25 | 51.960 | 61.955 | 21.718 | 27.136 |
| 0.5 | 53.072 | 63.736 | 25.190 | 32.123 |
| 1 | 54.768 | 65.350 | 29.598 | 38.458 |
| 2 | 57.603 | 67.961 | 35.073 | 43.194 |
| 5 | 62.909 | 74.836 | 43.778 | 52.991 |

All primary epsilon-zero entries indicate exact ambiguity. The full-pre
diagnostic also retains exact ambiguity. The comparison documents the value of
more historical geometry but must not be reinterpreted as a feasible population
certificate.

## 7. Static benchmarks and fairness

The role-matched outcome-only SCM is the fair comparison: 1975--1985 fitting
RMSE 0.2837 and 1986--1988 validation RMSE 1.0711. The full-pre outcome-only
SCM uses more direct information and has 1970--1988 RMSE 1.6564. The published
ADH-weight trajectory has full-pre RMSE 1.7576, but the original predictor-based
optimization was not rerun.

The completed outputs are consistent with the frozen donor universe, split,
candidate ordering, R grid, and separate q analyses. No output indicates
post-treatment treated-outcome tuning. The original fairness CSV is stale only
because it was not regenerated after R execution; its design checks remain
consistent with the scripts and output timestamps.

## 8. Observed treated outcomes

At q=1 and R=0, California's realized treated series lies below the plug-in
counterfactual interval in all 12 post-treatment years. At q=2's minimum radius
it lies inside three coordinate intervals and below nine. These comparisons
are descriptive. They do not convert realized outcomes into conditional means,
and they played no role in selecting q, kappa, h, R, donors, or certificates.

## 9. Inference classification

- **Tier 1:** available. The reported dynamic sets and fibers are plug-in
  identified sets under the supplied bridge and continuation assumptions.
- **Tier 2:** not run. No explicit maintained statistical allowances were
  supplied.
- **Tier 3:** unavailable. There is no defensible Gate-5 certificate bundle.

Specifically absent are a path-risk/anchor certificate, noise covariance and
restricted coercivity bounds, fitting-to-cutoff allowance, predictive-to-causal
transfer allowance, all-chord geometry epsilon, probability allocation, and a
joint q=2 level--velocity uncertainty region. No mixing/dependence input is
provided. The conic rational-dual certificates validate numerical upper
enclosures conditional on the plug-in optimization problems only.

## 10. Solver audit

No structural infeasibility is reported. The full requested coordinate endpoint
counts are present: 336 for each q, 96 for fibers, and 144 for deployment
moduli. Hence no failed coordinate endpoint was silently dropped.

- q=1 dynamic: 29 optimal, 297 close-to-optimal, 10 numerical-problem; all 336
  use rational-dual-residual certificates. Maximum recorded primal violation is
  0.111. Approximate primal values are not certified lower bounds because
  `primal_exact` is false.
- q=2 dynamic: 117 optimal, 191 close-to-optimal, 17 numerical-problem, 10
  multiplier-cone, and one maximum-iteration status. There are 323 rational
  dual certificates and 13 simplex-hull fallbacks. Some fallback rows have very
  large primal violations, up to 393.6; the certified hull endpoint, not the
  infeasible primal candidate, is the usable output.
- Evaluation fibers: 96/96 optimal with rational-dual certificates.
- Deployment modulus: 144/144 optimal with rational-dual certificates.

Path-functional support is less stable than coordinate support. q=1 total-path
support has two multiplier-cone and one numerical-problem statuses. q=2 total
support has four multiplier-cone and six numerical-problem statuses. These
fallback-heavy total-path values should remain supplementary.

## 11. Answers to the scientific questions

1. **Does q=1 produce an informative continuation set?** Yes, conditionally on
   the plug-in anchor and structural radius. Median width is 41.76 at R=0 and
   remains only 35.4% of the hull at the largest-grid median.
2. **Does q=2 produce an informative set?** At small R, yes, but the coordinate
   enclosures have material numerical degradation. At large R, the long-horizon
   set approaches the donor hull.
3. **Is q=2 wider because of velocity uncertainty?** Not in these outputs.
   Velocity sampling uncertainty was not propagated. q=2's widening is from
   second-order continuation slack and horizon accumulation; honest velocity
   uncertainty remains uncomputed.
4. **What is the main limitation?** Under the primary design, historical
   predictive ambiguity becomes substantial deployment ambiguity. q=1
   continuation uncertainty adds comparatively little; q=2 continuation
   uncertainty becomes important at moderate and large R. Across both, the
   absence of statistical certificates prevents honest inference.
5. **What does TRSC reveal that static SCM hides?** A single static trajectory
   masks 30--48-pack evaluation fibers and a median exact deployment modulus of
   50.62 packs.
6. **Does good fit imply identification?** No. Both TRSC fits have RMSE below
   \(7\times10^{-6}\), while counterfactual ambiguity remains economically large.
7. **Application or stress test?** It is a canonical application used most
   successfully as a partial-identification stress test.

## 12. Figures and tables

Recommended main-text set:

1. `figure_main_q1_counterfactual.png`: observed outcomes, fair static SCM, and
   the q=1 R=0 plug-in set. The caption must state that observed post-treatment
   outcomes are descriptive realizations.
2. `figure_main_dynamic_sensitivity.png`: q=1/q=2 sensitivity using the
   conservative monotone closure of certified endpoint enclosures and the
   q-specific multiplier c.
3. `figure_main_evaluation_deployment.png`: evaluation-fiber and deployment
   modulus diagnostics, with full-pre geometry explicitly labeled diagnostic.
4. `table_california_main_audited.csv`: compact quantitative sensitivity table.

Appendix:

- `figure1_static_benchmarks.png` for the published-ADH comparison;
- `figure2_trsc_pre_paths.png` for the interpolation result;
- `figure_geometry_spectrum.png` for rank diagnostics;
- raw `figure_q1_width_vs_R.png` and especially `figure_q2_width_vs_R.png` only
  with a warning that numerical enclosures need not be monotone. The raw q=2
  figure should not appear in the main text.

## 13. Final verdict

**CALIFORNIA IS PRIMARILY A PARTIAL-IDENTIFICATION STRESS TEST.**

Quantitatively, q=1's minimum-radius median width is 41.76 packs, q=2's is
35.09 before accounting for a fallback-driven maximum, the primary exact
deployment modulus has median 50.62, and no Tier-2 or Tier-3 statistical set is
available. The application is meaningful and publication-worthy if presented
as evidence about hidden ambiguity and deployment instability, not as a precise
new estimate of Proposition 99's effect.

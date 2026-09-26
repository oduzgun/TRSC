# California Proposition 99 TRSC results

## 1. Data and frozen design

The application uses the canonical annual cigarette-sales panel for 39 states
over 1970--2000. California is treated beginning in 1989 and the remaining 38
states form the donor pool. There are 19 untreated years and 12 deployment
years, with no missing cigarette-sales observations. The frozen role allocation
is calibration 1970--1974, fitting 1975--1985, and validation 1986--1988. The
equal-kernel local-linear validation weights targeted at 1988 are
\(-1/6,1/3,5/6\). No post-treatment California outcome entered tuning.

The original `run_manifest.json` and several status/placeholder files predate
the successful R run and incorrectly record TRSC as unexecuted. They are
preserved. The completed-run reconciliation is in
`output/run_manifest_completed_audit.json` and
`output/audit_reconciliation.csv`.

The terminal record also shows that the successful macOS run applied two
manual source fixes before execution: the labelled Stata state variable was
decoded with `haven::as_factor()`, and `trsc_boundary_readouts()` was added and
exported before reinstalling TRSC 0.1.0. The exact checksum of that executed
post-patch package tree was not captured. This limits byte-for-byte provenance
but does not alter the numerical interpretation below; the saved q=2 h=1
boundary agrees exactly with the adjacent 1985 and 1984 fitted coordinates.

## 2. Static SCM benchmarks

The primary fair static benchmark uses the same 1975--1985 fitting block and
38 donors as TRSC. Its fitting RMSE is 0.2837 and its 1986--1988 validation
RMSE is 1.0711. The descriptive full-pre outcome-only SCM uses 1970--1988 and
has full-pre RMSE 1.6564. Applying the published ADH weights gives full-pre
RMSE 1.7576; the original predictor-weight optimization was not re-estimated.
The predictor-based ADH design remains separate and predictor rows were not
treated as additional time observations.

## 3. TRSC pre-treatment selection

The structural orders were selected separately and were never compared to
choose a single causal order.

| q | selected candidate | kappa | h | validation score | fitting RMSE |
|---:|---|---:|---:|---:|---:|
| 1 | m019 | 30 | 1 | 2.4528 | \(6.47\times10^{-6}\) |
| 2 | m004 | 0.1 | 1 | 4.5415 | \(1.70\times10^{-6}\) |

For q=1, \(\kappa=30\) was the upper endpoint of the original grid. The
pre-treatment-only rule expanded the grid to 300, after which 30 remained the
interior minimizer. q=2 selected 0.1 on the original grid. No full-pre refit
was performed. Validation residuals for the selected q=1 readout are
-0.3551, -1.0925, and -1.5783 in 1986--1988; the corresponding q=2 residuals
are -0.5326, -1.3117, and -2.1951.

## 4. q=2 boundary state

The q=2 level and preceding readouts come from the same selected q=2 fit and
the same bandwidth \(h=1\). Their fitted endpoints are 1985 and 1984. In the
five primary metric coordinates,

\[
\widehat\beta_{-1}=(45.6207,47.0228,49.0360,49.6604,50.8848),
\]
\[
\widehat\beta_0=(45.1783,46.5730,48.5545,49.1497,50.3637),
\]
and
\[
\widehat v_0=(-0.4423,-0.4499,-0.4815,-0.5106,-0.5212),
\qquad \|\widehat v_0\|_2=1.0781.
\]

These are frozen predictive readouts. Because fitting ends in 1985, using them
as plug-in estimates of the 1988/1987 causal boundary requires an unquantified
fitting-to-cutoff transfer. No such transfer certificate is available.

## 5. Minimum radii and sensitivity grid

For q=1, persistence is analytically feasible:

\[
[R_{\min,1}^{L},R_{\min,1}^{U}]=[0,0].
\]

For q=2, the solver returned

\[
[R_{\min,2}^{L},R_{\min,2}^{U}]
=[0,1.48\times10^{-8}],
\]

with optimal status. Thus affine continuation is feasible to numerical
precision; California does not require a substantively positive minimum
curvature radius.

The scale factors are \(S_{1,H}=4.0734\times10^{-5}\) and
\(S_{2,H}=0.0028043\). Starting from the prespecified multipliers
\(c\in\{0,.25,.5,1,2,4,8,16\}\), the fixed doubling rule extended both grids
through \(c=1024\), because neither analysis reached donor-hull saturation.

## 6. Dynamic plug-in identified sets

At the minimum radius, q=1 widths range from 32.52 to 47.74 packs, with median
41.76. The median width is 33.7% of the unrestricted donor-hull width. Even at
the largest examined radius, \(R=0.04171\), the median width is 43.81 and the
largest hull fraction is 41.6%. The q=1 continuation restriction therefore
continues to tighten the donor hull, but most of its remaining width is already
present at \(R=0\) through deployment evaluation fibers.

At the q=2 minimum radius, the median width is 35.09. The reported maximum of
97.70 is not a sharp q=2 result: the 1993 upper endpoint used a certified
donor-hull fallback after the conic solver lost a usable dual. As curvature
slack grows, q=2 widens materially. At \(c=1024\), \(R=2.8716\), the median
width is 79.68 packs and the 2000 width is 106.73, or 92.4% of the donor hull.
The grid still did not formally saturate.

The exact identified sets are nested in R. The raw q=1 certified enclosures
have only three negligible numerical monotonicity deviations. The raw q=2
enclosures have 29 deviations, including large movements caused by donor-hull
fallbacks. The audited sensitivity figure therefore uses a conservative
monotone closure of the certified endpoint enclosures; the unmodified raw
outputs remain preserved.

## 7. Evaluation ambiguity

Using the primary 1970--1974 calibration geometry, no post-treatment date is
point identified. q=1 evaluation-fiber widths range from 32.52 to 47.74 packs
(median 41.75), and q=2 level-profile widths range from 30.03 to 42.46
(median 37.26). These widths are about one-third of the corresponding donor
hulls.

The full-pre geometry is diagnostic only. It reduces the median q=1 fiber to
2.66 and the median q=2 fiber to 6.83, but no date becomes point identified.
This comparison shows how much the result depends on the information assigned
to geometry calibration; it is not evidence that the full-pre metric is the
population geometry.

Rank deficiency establishes nonunique donor representations, not automatically
nonunique deployment values. Here nonunique deployment evaluation is established
separately by the positive optimized fiber widths.

## 8. Deployment amplification

Under the primary geometry, the exact deployment modulus \(\omega_x(0)\) ranges
from 41.10 to 59.51 packs, with median 50.62, approximately 41% of the donor
hull. Median modulus values rise to 54.77 at \(\epsilon=1\) and 62.91 at
\(\epsilon=5\). Thus relationships that are exactly or nearly equivalent on
the five calibration dates can imply materially different post-treatment
evaluations.

Under the full-pre diagnostic geometry, median \(\omega_x(0)\) is smaller but
still positive at 17.84. At \(\epsilon=5\), its median is 43.78. More historical
geometry narrows exact fibers, but deployment amplification remains visible.

## 9. Inference status

The empirical results are Tier 1 plug-in identified sets only. Tier 2 was not
run because no maintained numerical allowances were supplied. Tier 3 honest
Gate-5 inference is unavailable. Missing inputs include the path-risk envelope,
noise covariance and restricted coercivity bounds, fitting-to-cutoff and
predictive-to-causal transfer allowances, an all-simplex-chord geometry bound,
and valid probability allocation. A joint q=2 level--velocity statistical
certificate and dependence inputs are consequently absent as well.

The conic endpoints have numerical upper certificates conditional on their
mathematical programs. Those numerical certificates do not turn the plug-in
sets into statistical confidence regions.

## 10. Solver audit

Both OSQP fits report solved status. All 96 evaluation-fiber endpoints and all
144 deployment-modulus programs report optimal solutions with rational-dual
certificates.

The q=1 dynamic log contains all 336 requested endpoints: 29 optimal, 297 close
to optimal, and 10 numerical-problem statuses. Every endpoint nevertheless has
a rational-dual-residual upper certificate. The q=2 log also contains all 336
endpoints: 117 optimal, 191 close to optimal, 17 numerical-problem, 10
multiplier-cone, and one maximum-iteration status. It uses 323 rational-dual
certificates and 13 simplex-hull fallbacks. No structural infeasibility occurred
and no failed endpoint was dropped.

## 11. Static SCM comparison and observed outcomes

Static SCM produces a single precise-looking trajectory. TRSC shows that this
precision does not survive once predictively equivalent donor lifts and
assumption-indexed continuation are retained. The observed California series
lies below the q=1 minimum-radius plug-in set in all 12 post-treatment years;
it lies inside the q=2 minimum-radius outer interval in three years and below
it in nine. This is descriptive only: the realized treated outcome is neither
the untreated conditional mean nor a tuning target.

## 12. Scientific conclusion

California strongly demonstrates the separation between pre-treatment fit and
post-treatment identification. Both TRSC orders nearly interpolate the fitting
block, yet primary-geometry deployment fibers remain 30--48 packs wide and the
median exact deployment modulus is 50.62 packs. q=1 supplies stable and
meaningful tightening relative to the donor hull. q=2 is informative only at
small curvature radii and becomes much less informative with increasing R;
its current raw sensitivity output is also numerically less stable. Because no
Tier-2 or Tier-3 inference is available, the application cannot support claims
of honest sampling precision.

**Final verdict: CALIFORNIA IS PRIMARILY A PARTIAL-IDENTIFICATION STRESS TEST.**

It remains scientifically useful and canonical, but its strongest evidence is
about evaluation ambiguity, deployment amplification, and the fragility of
precision—not about a narrow estimated treatment effect.

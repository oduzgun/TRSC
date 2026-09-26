# California Proposition 99

We apply TRSC to the canonical California Proposition 99 cigarette-sales
panel. The data contain annual per-capita cigarette sales for California and
38 donor states from 1970 through 2000; treatment begins in 1989. We freeze the
empirical design before examining post-treatment outcomes. Donor-only geometry
is calibrated on 1970--1974, the predictive paths are fitted on 1975--1985,
and candidates are selected on 1986--1988 using local-linear weights targeted
at the 1988 cutoff. The static outcome-only benchmark uses the same donor pool
and fitting block. We analyze structural orders q=1 and q=2 separately rather
than selecting a causal continuation law from post-treatment performance.

For q=1, validation selects \(\kappa=30\) and terminal bandwidth \(h=1\).
Because 30 was the upper endpoint of the original penalty grid, the
pre-treatment-only rule expanded the grid to 300; 30 remained the interior
minimizer. For q=2, validation selects \(\kappa=0.1\) and \(h=1\) without grid
expansion. The fitting RMSEs are below \(7\times10^{-6}\) for both orders.
These nearly exact fits do not imply that the donor weights or post-treatment
evaluations are identified. They instead reflect the flexibility of the
rank-deficient short panel.

The q=2 boundary readouts use the same frozen fit and the same taper at adjacent
fitted endpoints. In the five calibration coordinates, the terminal velocity
has Euclidean norm 1.078. The minimum first-order continuation radius is
exactly zero. The certified q=2 minimum-radius interval is
\([0,1.48\times10^{-8}]\), so affine continuation is feasible to numerical
precision. These are plug-in boundary estimates based on fitting endpoints
1985 and 1984. The application does not supply the statistical transfer bound
needed to certify their relationship to the 1988 causal boundary.

At the minimum continuation radius, the q=1 coordinate bounds have median
width 41.76 packs per capita, or 33.7% of the unrestricted donor-hull width.
The q=1 widths change comparatively little over the prespecified sensitivity
grid: even after the fixed doubling rule extends the multiplier to 1024, the
median width is 43.81 and the largest horizon-specific hull fraction is 41.6%.
Thus q=1 supplies meaningful dynamic tightening, although its intervals are not
narrow in absolute terms.

The q=2 results are more radius-sensitive. At the minimum radius, the median
reported width is 35.09 packs. One maximum-radius statistic at the minimum
radius is driven by a donor-hull numerical fallback and should not be read as a
sharp support value. As the curvature allowance increases, q=2 expands much
more rapidly than q=1. At the largest reported radius, the median width is
79.68 packs and the year-2000 width is 106.73, equal to 92.4% of the donor
hull. The exact sets are nested in the radius, but several raw q=2 endpoint
enclosures are numerically nonmonotone because the solver falls back to valid
but loose donor-hull certificates. We therefore report a conservative monotone
closure in the main sensitivity figure and retain the raw solver output in the
replication archive.

Most of the q=1 width at zero continuation radius is evaluation rather than
continuation uncertainty. Holding the selected profile fixed, the primary
calibration geometry produces evaluation-fiber widths between 32.52 and 47.74
packs, with median 41.75. The corresponding q=2 level-profile median is 37.26.
No deployment date is point identified. A diagnostic calculation based on the
full 1970--1988 donor history reduces the median fibers to 2.66 for q=1 and
6.83 for q=2, but that geometry uses information outside the primary
calibration role and is not treated as population truth.

Deployment amplification is likewise substantial. Under the primary geometry,
the modulus at exact historical equivalence, \(\omega_x(0)\), ranges from
41.10 to 59.51 packs across post-treatment years and has median 50.62. At
predictive distance \(\epsilon=5\), the median modulus is 62.91. Hence two
relationships that are indistinguishable, or nearly indistinguishable, on the
calibration states can imply economically different evaluations at later donor
states. This is the central distinction from a static SCM trajectory: a single
well-fitting convex combination suppresses both the nonuniqueness of predictive
representations and their deployment consequences.

For reference, the role-matched static SCM has fitting RMSE 0.284 and
validation RMSE 1.071. The full-pre outcome-only SCM has 1970--1988 RMSE 1.656.
The realized California post-treatment series lies below the q=1 minimum-radius
plug-in interval in every year. We treat that comparison as descriptive only;
the realized treated outcome is not the untreated conditional mean and was not
used to choose q, the penalty, the readout bandwidth, the sensitivity radius,
or the donor pool.

The reported TRSC regions are plug-in identified sets, not confidence regions.
The application supplies valid numerical upper enclosures for the optimization
problems, including explicit donor-hull fallbacks when a useful conic dual is
unavailable. It does not supply a path-risk certificate, fitting-to-cutoff or
predictive-to-causal transfer bounds, an all-chord geometry certificate, or a
joint q=2 level--velocity uncertainty region. Tier-2 diagnostic outer sets and
Tier-3 honest Gate-5 confidence sets are therefore not reported.

California consequently functions as a partial-identification stress test. It
shows that excellent pre-treatment fit can coexist with large deployment
ambiguity, that q=1 continuation remains informative but not precise, and that
q=2 becomes weak at long horizons when curvature slack is allowed. The
application's contribution is not that TRSC produces a more favorable point
estimate than classical SCM. It is that the apparent precision of the
classical trajectory depends on suppressing evolving relationships,
nonunique evaluation fibers, and statistical uncertainty that the TRSC
framework makes explicit.

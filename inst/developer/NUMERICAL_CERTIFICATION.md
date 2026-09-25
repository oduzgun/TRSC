# Numerical input and optimization contract

## Mathematical program versus solver candidate

The conic backend maximizes c'z subject to Az=b, Gz+s=h, s in a product of
nonnegative and second-order cones. Variables are partitioned into unit simplex
blocks, with optional explicitly bounded scalar variables. ECOS is only a
candidate generator. Infeasible or inaccurate solver status is retained.

For any free y and cone-feasible d, write r=c-A'y-G'd. Weak duality gives

    c'z <= b'y + h'd + support_product_of_simplices_and_boxes(r).

The final term is a sum of block maxima plus positive-part box contributions.
GMP exact rational arithmetic evaluates this expression at the intended affine
coefficients. Nonnegative dual components are clipped upward to zero; each SOC
dual head is raised until its squared value is at least the exact tail squared
norm. No stationarity tolerance is treated as zero. Final conversion to a double
is checked to round outward. The unrestricted donor-hull support is another
valid bound and can be used when a solver dual is unusable.

These are certificates for coefficients supplied to the program (binary-double
inputs interpreted exactly). They are not interval certification of an unknown
population covariance, of a user's external floating-point preprocessing, or
of a causal assumption. Future-chord and numerical-anchor certificates must
cover any discrepancy between those intended population objects and the
recorded computational factor/centers.

## Feasibility, gaps and boundaries

Raw primal feasibility, raw dual objective, status and controls are retained.
An exactly verified feasible primal candidate gives a lower bound and a
certified gap. Otherwise the gap is NA; the approximate primal objective is
descriptive. No false zero gap is created from a solver convergence flag.

A negative certified support bound for the zero objective supplies a Farkas
infeasibility certificate. If the solver cannot certify emptiness, the return
can remain a conservative hull enclosure with an inconclusive solver status;
it must not be read as proof that a witness exists. In particular, a radius
between certified lower/upper Rmin bounds is numerically unresolved.

Profile equalities are encoded as two inequalities. Rank-deficient forward
factors and redundant equality rows do not require inverses or dropped weak
directions. Exact and minimum-radius faces require no Slater condition for
weak duality. No active-face or donor-subset enumeration is performed.

## Input responsibilities

| Input | Responsibility |
|---|---|
| Donor states, deployment x, sample dates | Observed, aligned and dimension checked |
| Empirical factor | Computed without positive-direction truncation |
| Structural q, R, auxiliary law | User-specified scientific restrictions |
| Anchor errors and all-chord epsilon | Maintained or independently certified numerical bounds |
| Risk envelope Bbar, fitting-event probability, gap and causal discrepancy | Supplied to the Markov helper; never estimated from unknown oracle quantities |
| Candidate events, probabilities and covered IDs | User-certified simultaneous scope; union accounting is performed by software |
| Mixing error | Supplied valid coupling allowance; no mixing estimation |
| W | User choice, checked positive definite; numerical factorization residue is explicitly enlarged |
| Solver dual and primal candidates | Computed; certificate checker validates/corrects independently of status |
| Treated mean functional | Known by assumption or separately certified; not a realized outcome |

Sum-one roundoff in a nonnegative representative is normalized using exact
rational proportions. Negative weights are rejected. Bounds must refer to the
recorded center; any change from a fitted readout needs a numerical allowance.
Two numerical minimizers may give different confidence endpoints; honest
coverage requires the supplied risk envelope to hold uniformly over admissible
approximate minimizers. A certificate constructor cannot verify that theorem.

## Performance

Dense exact rational matrix arithmetic is deliberately conservative and can
be costly at large J, H or calibration row count N. Complexity is polynomial
in program dimensions; bit growth and repeated directional solves can dominate.
There is no unverified low-rank compression or parallel subset search.
Future optimized backends must preserve this certification contract.

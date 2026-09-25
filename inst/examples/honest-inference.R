library(TRSC)

g <- trsc_geometry(matrix(c(0, 1), 1, dimnames = list(NULL, c("A", "B"))))
x <- rbind(c(A = 0, B = 1), c(A = 0, B = 1))

# Illustrative maintained bounds, not estimates or universal default constants.
# Each component is in the joint W norm for a q=2 boundary certificate.
anchor <- trsc_anchor_certificate(g, c(A = .5, B = .5), q = 2,
  previous = c(A = .6, B = .4), W = matrix(c(2, .5, .5, 1), 2),
  components = c(estimation = .05, fitting_gap = .02, causal_gap = 0, numerical = 0),
  failure_probability = .02, source = "Toy maintained joint boundary bounds")
geometry <- trsc_geometry_certificate(g, epsilon = .01,
  failure_probability = .01, source = "Toy maintained ALL-chord bound")

outer <- trsc_confidence_set(x, anchor, geometry, q = 2, R = .1)
print(outer)
print(trsc_bounds(outer))
print(trsc_anchor_propagation(anchor, 1:5))
support <- trsc_confidence_support(outer, c(1, 1))
stopifnot(support$certificate$valid,
          abs(outer$energy_allowance - 2 * sqrt(2) * .01) < 1e-12)
print(support)
print(trsc_diagnostics(support))

# 3 is a supplied conditional-MEAN functional, not an observed realization.
print(trsc_effect_interval(outer, c(1, 1), treated_mean = 3))

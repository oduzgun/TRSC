library(TRSC)

# Donors A and B agree predictively but differ at deployment.
g <- trsc_geometry(matrix(c(0, 0, 1), 1,
                         dimnames = list(NULL, c("A", "B", "C"))))
x <- c(A = 0, B = 2, C = 1)
fiber <- trsc_evaluation_set(g, x, representative = c(A = 1, B = 0, C = 0))
stopifnot(abs(trsc_evaluation_width(fiber) - 2) < 1e-6)

# R=0 fixes the profile, not its post-treatment evaluations.
z <- trsc_counterfactual_set(g, rbind(x, x), c(A = 1, B = 0, C = 0),
                            q = 1, R = 0, type = "population")
print(trsc_bounds(z))
print(trsc_support(z, c(1, -1)))

# A scientifically justified, unaffected measurement can remove the fiber.
za <- trsc_auxiliary_refine(z, maps = list(matrix(x, 1), matrix(x, 1)),
  values = list(.5, .5), source = "Toy exact untreated shared-composition measurement")
stopifnot(all(abs(trsc_bounds(za)$upper - .5) < 1e-6))

line <- trsc_geometry(matrix(c(0, 1), 1, dimnames = list(NULL, c("A", "B"))))
minimum <- trsc_min_radius(line, c(A = 0, B = 1), q = 2, H = 2,
                           previous = c(A = 1, B = 0))
stopifnot(abs(minimum$lower - 1) < 1e-6, abs(minimum$upper - 1) < 1e-6)
print(minimum[c("lower", "upper", "estimate", "status")])

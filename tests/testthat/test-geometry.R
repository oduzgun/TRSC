test_that("duplicate donors are allowed and profile-equivalent", {
  X <- cbind(
    A = c(1, 2, 3, 4, 5),
    B = c(2, 1, 0, 1, 2),
    A_duplicate = c(1, 2, 3, 4, 5)
  )
  g <- trsc_geometry(X)

  expect_s3_class(g, "trsc_geometry")
  expect_equal(g$Gamma, crossprod(X) / nrow(X), tolerance = 1e-12)
  expect_false(g$diagnostics$representation_unique_empirical)
  expect_gte(g$diagnostics$empirical_null_dimension, 1L)

  w <- c(A = 1, B = 0, A_duplicate = 0)
  v <- c(A = 0, B = 0, A_duplicate = 1)
  expect_equal(trsc_profile_distance(g, w, v), 0, tolerance = 1e-10)
  expect_true(trsc_profile_equivalent(g, w, v))
  expect_equal(
    trsc_profile_coordinates(g, w),
    trsc_profile_coordinates(g, v),
    tolerance = 1e-10
  )
})

test_that("donor permutation preserves predictive geometry", {
  X <- cbind(
    A = c(1, 2, 4, 8),
    B = c(2, 0, 3, 5),
    C = c(1, 3, 2, 7)
  )
  g1 <- trsc_geometry(X)
  perm <- c("C", "A", "B")
  g2 <- trsc_geometry(X[, perm, drop = FALSE])

  w <- c(A = 0.2, B = 0.3, C = 0.5)
  v <- c(A = 0.5, B = 0.4, C = 0.1)

  expect_equal(
    trsc_profile_distance(g1, w, v),
    trsc_profile_distance(g2, w[perm], v[perm]),
    tolerance = 1e-12
  )
})

test_that("singular geometry is a valid state, not an error", {
  X <- cbind(
    A = 1:6,
    B = 2 * (1:6),
    C = 3 * (1:6)
  )
  g <- expect_no_error(trsc_geometry(X))
  expect_lt(g$diagnostics$empirical_rank, ncol(X))
  expect_gt(g$diagnostics$empirical_null_dimension, 0L)
})

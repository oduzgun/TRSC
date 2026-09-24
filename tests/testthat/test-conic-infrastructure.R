test_that("rank diagnostics never truncate a positive predictive direction", {
  X <- cbind(A = c(0, 0), B = c(1e-10, -1e-10))
  g <- trsc_geometry(X, tolerance = 1)
  expect_equal(sum((g$B %*% c(-1, 1))^2), 1e-20, tolerance = 1e-30)
  expect_equal(trsc_profile_distance(g, c(1, 0), c(0, 1)), 1e-10,
               tolerance = 1e-20)
})

test_that("exact residual correction certifies an imperfect dual candidate", {
  p <- .trsc_conic_problem(2L, list(1:2))
  p <- .trsc_add_equality(p, matrix(c(1, 1), 1), 1)
  p <- .trsc_add_soc(p, matrix(c(0, 1), 1), 0, 0.4)
  ans <- .trsc_solve_socp(p, c(0, 1))
  expect_true(ans$certificate$valid)
  expect_gte(ans$upper, 0.4)
  expect_lt(ans$upper, 0.4001)
  bad <- .trsc_certify_upper(p, c(0, 1), y = 0.001,
                           z = c(0, 0, 0.1, 0.7))
  expect_gte(bad$upper, 0.4)
  expect_true(bad$valid)
})

test_that("a failed solver cannot turn a primal guess into an upper endpoint", {
  p <- .trsc_conic_problem(2L, list(1:2))
  p <- .trsc_add_equality(p, matrix(c(1, 1), 1), 1)
  ans <- .trsc_solve_socp(p, c(-2, 3), control = list(maxit = 1L))
  expect_true(ans$certificate$valid)
  expect_gte(ans$upper, 3)
  expect_true(is.list(ans$raw))
})

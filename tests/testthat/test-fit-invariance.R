make_trsc_panel <- function(duplicate = FALSE, order = NULL) {
  periods <- 1:7
  A <- c(1, 2, 4, 3, 6, 5, 8)
  B <- c(7, 5, 4, 6, 2, 3, 1)
  target <- 0.4 * A + 0.6 * B

  units <- if (duplicate) c("Target", "A", "B", "A2") else c("Target", "A", "B")
  if (!is.null(order)) units <- c("Target", order)
  panel <- expand.grid(unit = units, period = periods, stringsAsFactors = FALSE)

  lookup <- list(Target = target, A = A, B = B, A2 = A)
  panel$outcome <- unlist(lapply(periods, function(t) {
    vapply(units, function(u) lookup[[u]][t], numeric(1))
  }), use.names = FALSE)
  panel
}

test_that("duplicate donor representations leave fitted reference predictions unchanged", {
  base_panel <- make_trsc_panel(FALSE)
  dup_panel <- make_trsc_panel(TRUE)

  p1 <- trsc_data(base_panel, "unit", "period", "outcome", "Target",
                  c("A", "B"), 8, 1:7, "outcome")
  p2 <- trsc_data(dup_panel, "unit", "period", "outcome", "Target",
                  c("A", "B", "A2"), 8, 1:7, "outcome")

  ref <- cbind(
    A = c(-2, -1, 0, 1, 2, 3),
    B = c(3, 1, 2, -1, 0, 4)
  )
  ref_dup <- cbind(A = ref[, "A"], B = ref[, "B"], A2 = ref[, "A"])

  f1 <- trsc_fit(p1, trsc_geometry(ref), kappa = 3, q = 1)
  f2 <- trsc_fit(p2, trsc_geometry(ref_dup), kappa = 3, q = 1)

  eval1 <- ref %*% stats::weights(f1)
  eval2 <- ref_dup %*% stats::weights(f2)

  expect_equal(eval1, eval2, tolerance = 2e-5)
  expect_equal(stats::fitted(f1), stats::fitted(f2), tolerance = 2e-5)
  expect_gt(f2$diagnostics$simplex_null_dimension_empirical, 0L)
  expect_false(f2$diagnostics$representation_unique_under_empirical_geometry)
})

test_that("donor permutation changes coordinates but not the fitted predictive relationship", {
  panel1 <- make_trsc_panel(FALSE)
  panel2 <- make_trsc_panel(FALSE, order = c("B", "A"))

  p1 <- trsc_data(panel1, "unit", "period", "outcome", "Target",
                  c("A", "B"), 8, 1:7, "outcome")
  p2 <- trsc_data(panel2, "unit", "period", "outcome", "Target",
                  c("B", "A"), 8, 1:7, "outcome")

  ref <- cbind(
    A = c(-2, 0, 1, 3, 4),
    B = c(4, 1, -1, 0, 2)
  )
  ref_perm <- ref[, c("B", "A"), drop = FALSE]

  f1 <- trsc_fit(p1, trsc_geometry(ref), kappa = 2, q = 1)
  f2 <- trsc_fit(p2, trsc_geometry(ref_perm), kappa = 2, q = 1)

  eval1 <- ref %*% stats::weights(f1)
  eval2 <- ref_perm %*% stats::weights(f2)
  expect_equal(eval1, eval2, tolerance = 2e-5)
  expect_equal(stats::fitted(f1), stats::fitted(f2), tolerance = 2e-5)
})

test_that("singular geometry is accepted by the fitting program", {
  panel <- make_trsc_panel(TRUE)
  prepared <- trsc_data(panel, "unit", "period", "outcome", "Target",
                        c("A", "B", "A2"), 8, 1:7, "outcome")
  ref <- cbind(
    A = c(1, 2, 3, 4, 5),
    B = c(5, 3, 4, 1, 2),
    A2 = c(1, 2, 3, 4, 5)
  )
  geometry <- trsc_geometry(ref)

  expect_lt(geometry$diagnostics$empirical_rank, 3L)
  fit <- expect_no_error(trsc_fit(prepared, geometry, kappa = 5, q = 1))
  expect_s3_class(fit, "trsc_fit")
  expect_true(tolower(fit$solver$status) %in% c("solved", "solved inaccurate"))
})

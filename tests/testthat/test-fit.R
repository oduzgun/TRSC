test_that("trsc_fit recovers an exactly representable static relationship", {
  panel <- expand.grid(
    unit = c("Target", "A", "B"),
    period = 1:8,
    stringsAsFactors = FALSE
  )
  A <- seq(1, 8)
  B <- seq(10, 24, length.out = 8)
  target <- 0.3 * A + 0.7 * B
  values <- c(rbind(target, A, B))
  panel$outcome <- values

  prepared <- trsc_data(
    panel, "unit", "period", "outcome", "Target", c("A", "B"),
    intervention_time = 9, pre_period = 1:8, predictors = "outcome"
  )

  calibration <- cbind(
    A = c(-2, -1, 0, 1, 2, 3),
    B = c(3, 1, -1, 2, 4, 0)
  )
  geometry <- trsc_geometry(calibration)
  fit <- trsc_fit(prepared, geometry, kappa = 10, q = 1)

  expect_s3_class(fit, "trsc_fit")
  expect_equal(unname(stats::fitted(fit)), target, tolerance = 2e-5)
  expect_equal(unname(colSums(stats::weights(fit))), rep(1, 8), tolerance = 2e-6)
  expect_true(all(stats::weights(fit) >= -2e-6))
  expect_equal(fit$objective$fit_loss, 0, tolerance = 1e-8)
  expect_equal(fit$objective$temporal_roughness, 0, tolerance = 1e-8)
  expect_identical(fit$specification$reference_law_certified, FALSE)
})

test_that("q = 2 leaves an affine predictive relationship unpenalized", {
  periods <- 1:9
  donor_a <- rep(0, length(periods))
  donor_b <- rep(10, length(periods))
  w_b <- seq(0.1, 0.9, length.out = length(periods))
  target <- 10 * w_b

  panel <- expand.grid(
    unit = c("Target", "A", "B"),
    period = periods,
    stringsAsFactors = FALSE
  )
  panel$outcome <- c(rbind(target, donor_a, donor_b))

  prepared <- trsc_data(
    panel, "unit", "period", "outcome", "Target", c("A", "B"),
    intervention_time = 10, pre_period = periods, predictors = "outcome"
  )
  geometry <- trsc_geometry(cbind(
    A = c(0, 0, 0, 0),
    B = c(8, 10, 12, 14)
  ))

  fit <- trsc_fit(prepared, geometry, kappa = 1000, q = 2)
  expect_equal(unname(stats::fitted(fit)), target, tolerance = 5e-5)
  expect_equal(fit$objective$temporal_roughness, 0, tolerance = 5e-8)
})

test_that("fit requires an explicitly compatible supplied geometry", {
  panel <- expand.grid(
    unit = c("Target", "A", "B"),
    period = 1:5,
    stringsAsFactors = FALSE
  )
  panel$outcome <- seq_len(nrow(panel))
  prepared <- trsc_data(
    panel, "unit", "period", "outcome", "Target", c("A", "B"),
    6, 1:5, "outcome"
  )
  wrong <- trsc_geometry(cbind(A = 1:4, C = 4:1))
  expect_error(
    trsc_fit(prepared, wrong, kappa = 1),
    "must match the fitted donor set exactly"
  )
  expect_error(trsc_fit(prepared, trsc_geometry(cbind(A=1:4,B=4:1)),
                        kappa = -1), "nonnegative")
  expect_error(trsc_fit(prepared, trsc_geometry(cbind(A=1:4,B=4:1)),
                        kappa = 1, q = 3), "exactly 1 or 2")
})


test_that("predictive_path separates the fitted relationship from metric identification claims", {
  panel <- expand.grid(
    unit = c("Target", "A", "B"),
    period = 1:6,
    stringsAsFactors = FALSE
  )
  A <- 1:6
  B <- c(6, 4, 5, 2, 3, 1)
  target <- 0.25 * A + 0.75 * B
  panel$outcome <- c(rbind(target, A, B))

  prepared <- trsc_data(
    panel, "unit", "period", "outcome", "Target", c("A", "B"),
    7, 1:6, "outcome"
  )
  geometry <- trsc_geometry(cbind(
    A = c(-1, 0, 1, 2),
    B = c(2, -1, 3, 0)
  ))
  fit <- trsc_fit(prepared, geometry, kappa = 2, q = 1)

  path <- predictive_path(fit)
  expect_s3_class(path, "trsc_predictive_path")
  expect_identical(dim(path$representative), dim(stats::weights(fit)))
  expect_equal(as.numeric(path$representative),
               as.numeric(stats::weights(fit)), tolerance = 0)
  expect_equal(path$metric_coordinates, fit$metric_coordinates)
  expect_identical(
    path$diagnostics$representation_unique_under_empirical_geometry,
    fit$diagnostics$representation_unique_under_empirical_geometry
  )

  diagnostics <- trsc_diagnostics(fit)
  expect_s3_class(diagnostics, "trsc_diagnostics")
  expect_equal(diagnostics$fit$sse, fit$objective$fit_loss)
  expect_false(diagnostics$certification$population_coercivity)
  expect_false(diagnostics$certification$calibration_complexity)
  expect_false(diagnostics$certification$causal_identification)
})


test_that("matching-predictor values do not enter the TRSC outcome-regression objective", {
  periods <- 1:6
  A <- c(1, 3, 2, 5, 4, 7)
  B <- c(7, 5, 6, 2, 3, 1)
  target <- 0.4 * A + 0.6 * B

  panel <- expand.grid(
    unit = c("Target", "A", "B"),
    period = periods,
    stringsAsFactors = FALSE
  )
  panel$outcome <- c(rbind(target, A, B))
  panel$aux <- seq_len(nrow(panel))
  panel2 <- panel
  panel2$aux <- panel2$aux^3 + 17

  p1 <- trsc_data(
    panel, "unit", "period", "outcome", "Target", c("A", "B"),
    7, periods, c("outcome", "aux")
  )
  p2 <- trsc_data(
    panel2, "unit", "period", "outcome", "Target", c("A", "B"),
    7, periods, c("outcome", "aux")
  )
  geometry <- trsc_geometry(cbind(
    A = c(-2, -1, 1, 3, 4),
    B = c(4, 2, -1, 0, 3)
  ))

  f1 <- trsc_fit(p1, geometry, kappa = 4, q = 1)
  f2 <- trsc_fit(p2, geometry, kappa = 4, q = 1)
  expect_equal(stats::fitted(f1), stats::fitted(f2), tolerance = 1e-7)
  expect_equal(f1$objective$total, f2$objective$total, tolerance = 1e-7)
})

test_that("reported TRSC objective equals the direct Gamma-weighted objective", {
  periods <- 1:7
  A <- c(1, 2, 4, 3, 6, 5, 8)
  B <- c(8, 7, 5, 6, 2, 3, 1)
  target <- c(5.2, 5.1, 4.7, 4.9, 3.0, 3.2, 2.1)

  panel <- expand.grid(
    unit = c("Target", "A", "B"),
    period = periods,
    stringsAsFactors = FALSE
  )
  panel$outcome <- c(rbind(target, A, B))
  prepared <- trsc_data(
    panel, "unit", "period", "outcome", "Target", c("A", "B"),
    8, periods, "outcome"
  )
  geometry <- trsc_geometry(cbind(
    A = c(-1, 0, 1, 2, 4),
    B = c(3, 1, -2, 0, 2)
  ))

  for (order in c(1L, 2L)) {
    fit <- trsc_fit(prepared, geometry, kappa = 2.5, q = order)
    W <- unclass(stats::weights(fit))
    U <- t(W)
    delta <- if (order == 1L) {
      U[-1, , drop = FALSE] - U[-nrow(U), , drop = FALSE]
    } else {
      U[-c(1, 2), , drop = FALSE] -
        2 * U[-c(1, nrow(U)), , drop = FALSE] +
        U[-c(nrow(U) - 1, nrow(U)), , drop = FALSE]
    }
    direct_roughness <- sum((delta %*% geometry$Gamma) * delta)
    direct_loss <- sum((prepared$blocks$Y1 - rowSums(prepared$blocks$Y0 * U))^2)

    expect_equal(fit$objective$temporal_roughness, direct_roughness,
                 tolerance = 1e-7)
    expect_equal(fit$objective$fit_loss, direct_loss, tolerance = 1e-7)
    expect_equal(fit$objective$total,
                 direct_loss + 2.5 * direct_roughness,
                 tolerance = 1e-7)
  }
})


test_that("roughness decreases and residual loss increases with kappa", {
  periods <- 1:8
  donor_a <- rep(0, length(periods))
  donor_b <- rep(1, length(periods))
  target <- c(0.10, 0.90, 0.15, 0.85, 0.20, 0.80, 0.25, 0.75)

  panel <- expand.grid(
    unit = c("Target", "A", "B"),
    period = periods,
    stringsAsFactors = FALSE
  )
  panel$outcome <- c(rbind(target, donor_a, donor_b))

  prepared <- trsc_data(
    panel, "unit", "period", "outcome", "Target", c("A", "B"),
    intervention_time = 9, pre_period = periods, predictors = "outcome"
  )
  geometry <- trsc_geometry(cbind(
    A = c(0, 0, 0, 0, 0),
    B = c(0.5, 1, 1.5, 2, 2.5)
  ))

  for (order in c(1L, 2L)) {
    low <- trsc_fit(prepared, geometry, kappa = 0.01, q = order)
    high <- trsc_fit(prepared, geometry, kappa = 50, q = order)

    expect_lte(
      high$objective$temporal_roughness,
      low$objective$temporal_roughness + 1e-7
    )
    expect_gte(
      high$objective$fit_loss + 1e-7,
      low$objective$fit_loss
    )
  }
})

test_that("kappa zero recovers period-specific feasible outcome fits", {
  periods <- 1:6
  donor_a <- rep(0, length(periods))
  donor_b <- rep(1, length(periods))
  target <- c(0.15, 0.80, 0.30, 0.65, 0.45, 0.55)

  panel <- expand.grid(
    unit = c("Target", "A", "B"),
    period = periods,
    stringsAsFactors = FALSE
  )
  panel$outcome <- c(rbind(target, donor_a, donor_b))

  prepared <- trsc_data(
    panel, "unit", "period", "outcome", "Target", c("A", "B"),
    intervention_time = 7, pre_period = periods, predictors = "outcome"
  )
  geometry <- trsc_geometry(cbind(
    A = c(0, 0, 0, 0),
    B = c(0.5, 1, 1.5, 2)
  ))

  fit <- trsc_fit(prepared, geometry, kappa = 0, q = 1)
  expect_equal(stats::fitted(fit), target, tolerance = 2e-6)
  expect_equal(fit$objective$fit_loss, 0, tolerance = 1e-10)
})

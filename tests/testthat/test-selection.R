make_selection_fixture <- function() {
  periods <- 1:10
  A <- rep(0, length(periods))
  B <- rep(10, length(periods))
  wB <- c(0.2, 0.2, 0.2, 0.8, rep(0.8, 6))
  target <- 10 * wB

  panel <- expand.grid(
    unit = c("Target", "A", "B"),
    period = periods,
    stringsAsFactors = FALSE
  )
  panel$outcome <- c(rbind(target, A, B))

  prepared <- trsc_data(
    panel, "unit", "period", "outcome", "Target", c("A", "B"),
    intervention_time = 11, pre_period = periods, predictors = "outcome"
  )
  geometry <- trsc_geometry(cbind(
    A = c(0, 0, 0, 0),
    B = c(8, 10, 12, 14)
  ))
  fit <- trsc_fit(
    prepared, geometry, kappa = 0, q = 1, periods = 1:4
  )
  list(prepared = prepared, geometry = geometry, fit = fit)
}

test_that("terminal taper matches the paper and remains a convex simplex readout", {
  fixture <- make_selection_fixture()
  readout <- trsc_readout(fixture$fit, h = 4)

  expect_s3_class(readout, "trsc_readout")
  expect_equal(
    unname(readout$taper_weights),
    c(0.4, 0.3, 0.2, 0.1),
    tolerance = 1e-14
  )
  expect_equal(sum(readout$taper_weights), 1, tolerance = 1e-14)
  expect_equal(sum(readout$representative), 1, tolerance = 2e-6)
  expect_true(all(readout$representative >= -2e-6))
  expect_true(readout$frozen)

  endpoint <- trsc_readout(fixture$fit, h = 1)
  expect_equal(
    as.numeric(endpoint$representative),
    as.numeric(stats::weights(fixture$fit)[, 4]),
    tolerance = 2e-6
  )
})

test_that("boundary readouts use one frozen fit and the same shifted taper", {
  fixture <- make_selection_fixture()
  boundary <- trsc_boundary_readouts(fixture$fit, h = 2)
  W <- stats::weights(fixture$fit)
  taper <- c(2 / 3, 1 / 3)

  expect_s3_class(boundary, "trsc_boundary_readouts")
  expect_s3_class(boundary$level, "trsc_readout")
  expect_s3_class(boundary$previous, "trsc_readout")
  expect_identical(boundary$level, trsc_readout(fixture$fit, h = 2))
  expect_equal(
    unname(boundary$level$representative),
    unname(drop(W[, c(4, 3), drop = FALSE] %*% taper)),
    tolerance = 2e-6
  )
  expect_equal(
    unname(boundary$previous$representative),
    unname(drop(W[, c(3, 2), drop = FALSE] %*% taper)),
    tolerance = 2e-6
  )
  expect_identical(boundary$level$cutoff, fixture$fit$periods[4])
  expect_identical(boundary$previous$cutoff, fixture$fit$periods[3])
  expect_equal(
    boundary$velocity_metric_coordinates,
    boundary$level$metric_coordinates - boundary$previous$metric_coordinates,
    tolerance = 1e-14
  )
  deployment <- matrix(
    c(A = 0, B = 10), nrow = 1,
    dimnames = list(NULL, c("A", "B"))
  )
  q2_set <- trsc_counterfactual_set(
    fixture$geometry, deployment,
    anchor = boundary$level,
    q = 2, R = 0,
    previous = boundary$previous
  )
  expect_s3_class(q2_set, "trsc_identified_set")
  expect_equal(q2_set$anchors[[1]], boundary$previous$representative,
               tolerance = 2e-6)
  expect_equal(q2_set$anchors[[2]], boundary$level$representative,
               tolerance = 2e-6)
  expect_true(boundary$frozen)
  expect_error(
    trsc_boundary_readouts(fixture$fit, h = 4),
    "one additional fitted date"
  )
})

test_that("default readout library uses dyadic bandwidths plus full length", {
  fixture <- make_selection_fixture()
  library <- trsc_readout_library(fixture$fit)

  expect_s3_class(library, "trsc_readout_library")
  expect_identical(names(library), c("h1", "h2", "h4"))
  expect_identical(
    vapply(library, function(x) x$bandwidth, integer(1)),
    c(h1 = 1L, h2 = 2L, h4 = 4L)
  )
})

test_that("cutoff-targeted selection enforces the local-linear moment conditions", {
  fixture <- make_selection_fixture()
  candidates <- trsc_readout_library(fixture$fit, bandwidths = c(1, 4))

  selection <- trsc_select(
    fixture$prepared,
    candidates,
    validation_periods = 5:9,
    target_period = 10
  )

  expect_s3_class(selection, "trsc_selection")
  expect_equal(sum(selection$validation_weights), 1, tolerance = 1e-12)
  expect_equal(
    sum(selection$validation_weights * (5:9 - 10)),
    0,
    tolerance = 1e-12
  )
  expect_true(selection$diagnostics$signed_weights)
  expect_true(selection$diagnostics$candidates_frozen)
  expect_false(selection$diagnostics$refit_performed)
  expect_false(selection$diagnostics$theorem_assumptions_certified)
})

test_that("frozen validation selects the correct endpoint readout in a noiseless design", {
  fixture <- make_selection_fixture()
  candidates <- trsc_readout_library(fixture$fit, bandwidths = c(1, 4))

  selection <- trsc_select(
    fixture$prepared,
    candidates,
    validation_periods = 5:9
  )

  expect_identical(selection$selected_name, "h1")
  expect_equal(selection$scores[["h1"]], 0, tolerance = 1e-8)
  expect_gt(selection$scores[["h4"]], 0)
  expect_identical(selection$selected, candidates[["h1"]])
})

test_that("selection rejects candidates that were fitted on validation outcomes", {
  fixture <- make_selection_fixture()
  contaminated_fit <- trsc_fit(
    fixture$prepared, fixture$geometry,
    kappa = 0, q = 1, periods = 1:6
  )
  contaminated <- list(h1 = trsc_readout(contaminated_fit, 1))

  expect_error(
    trsc_select(
      fixture$prepared,
      contaminated,
      validation_periods = 5:9
    ),
    "fitting cutoff must precede all validation periods"
  )
})

test_that("selection requires the last untreated date as the theorem-aligned target", {
  fixture <- make_selection_fixture()
  candidates <- list(h1 = trsc_readout(fixture$fit, 1))

  expect_error(
    trsc_select(
      fixture$prepared,
      candidates,
      validation_periods = 5:9,
      target_period = 9
    ),
    "must equal the last declared pre-treatment date"
  )
})


test_that("selection requires all candidates to share the same frozen fitting block", {
  fixture <- make_selection_fixture()
  fit_short <- fixture$fit
  fit_long <- trsc_fit(
    fixture$prepared, fixture$geometry,
    kappa = 0, q = 1, periods = 1:5
  )
  candidates <- list(
    short = trsc_readout(fit_short, 1),
    long = trsc_readout(fit_long, 1)
  )

  expect_error(
    trsc_select(
      fixture$prepared,
      candidates,
      validation_periods = 6:9
    ),
    "same fitting block and cutoff"
  )
})


test_that("exact validation ties use the frozen candidate order deterministically", {
  fixture <- make_selection_fixture()
  first <- trsc_readout(fixture$fit, 1)
  candidates <- list(first = first, second = first)

  selection <- trsc_select(
    fixture$prepared,
    candidates,
    validation_periods = 5:9
  )

  expect_identical(selection$selected_name, "first")
  expect_identical(selection$diagnostics$exact_tie_count, 2L)
  expect_identical(
    selection$diagnostics$exact_tie_names,
    c("first", "second")
  )
  expect_identical(
    selection$diagnostics$tie_rule,
    "first exact minimum in candidate order"
  )
})

make_workflow_fixture <- function() {
  periods <- 1:15
  A <- rep(0, length(periods))
  B <- rep(10, length(periods))
  weight_b <- c(
    rep(0.2, 4),
    0.2,
    0.2, 0.2, 0.2, 0.8,
    0.8,
    rep(0.8, 5)
  )
  target <- 10 * weight_b

  panel <- expand.grid(
    unit = c("Target", "A", "B"),
    period = periods,
    stringsAsFactors = FALSE
  )
  panel$outcome <- c(rbind(target, A, B))

  prepared <- trsc_data(
    panel,
    "unit", "period", "outcome",
    "Target", c("A", "B"),
    intervention_time = 16,
    pre_period = periods,
    predictors = "outcome"
  )

  plan <- trsc_sample_plan(
    prepared,
    calibration_periods = 1:4,
    fitting_periods = 6:9,
    validation_periods = 11:14,
    calibration_thin = 2,
    validation_thin = 2,
    target_period = 15,
    min_gap0 = 2,
    min_gap1 = 2
  )

  list(panel = panel, prepared = prepared, plan = plan)
}

test_that("sample plan implements deterministic thinning, intact fitting and gaps", {
  fixture <- make_workflow_fixture()
  plan <- fixture$plan

  expect_s3_class(plan, "trsc_sample_plan")
  expect_identical(plan$calibration$retained_periods, c(1L, 3L))
  expect_identical(plan$fitting$periods, 6:9)
  expect_identical(plan$validation$retained_periods, c(11L, 13L))
  expect_identical(plan$gaps$g0_steps, 3L)
  expect_identical(plan$gaps$g1_steps, 2L)
  expect_true(plan$fitting$intact)
  expect_true(plan$diagnostics$calibration_count_matches_formula)
  expect_true(plan$diagnostics$validation_count_matches_formula)
  expect_false(plan$diagnostics$theorem_assumptions_certified)
})

test_that("sample plan rejects role overlap, noncontiguous fitting and inadequate gaps", {
  fixture <- make_workflow_fixture()
  p <- fixture$prepared

  expect_error(
    trsc_sample_plan(
      p, 1:4, 4:8, 10:14,
      target_period = 15
    ),
    "calibration segment must end before"
  )
  expect_error(
    trsc_sample_plan(
      p, 1:3, c(5, 6, 8), 10:14,
      target_period = 15
    ),
    "one contiguous raw date segment"
  )
  expect_error(
    trsc_sample_plan(
      p, 1:4, 5:8, 9:14,
      target_period = 15,
      min_gap0 = 2
    ),
    "smaller than min_gap0"
  )
})

test_that("calibration geometry uses retained donor states and no target outcomes", {
  fixture <- make_workflow_fixture()
  geometry <- trsc_calibration_geometry(fixture$prepared, fixture$plan)

  expected <- fixture$prepared$blocks$Y0[
    fixture$plan$calibration$retained_index, , drop = FALSE
  ]
  expect_equal(geometry$Gamma, crossprod(expected) / nrow(expected),
               tolerance = 1e-14)
  expect_identical(geometry$provenance$retained_periods, c(1L, 3L))
  expect_false(geometry$provenance$target_outcomes_used)
  expect_false(geometry$provenance$fitting_outcomes_used)
  expect_false(geometry$provenance$validation_outcomes_used)
  expect_false(geometry$provenance$population_reference_law_certified)
})

test_that("candidate library freezes each q-kappa fit once and expands readouts", {
  fixture <- make_workflow_fixture()
  library <- trsc_candidate_library(
    fixture$prepared,
    fixture$plan,
    q = c(1, 2),
    kappa = c(0.05, 2),
    bandwidths = c(1, 2, 4)
  )

  expect_s3_class(library, "trsc_candidate_library")
  expect_length(attr(library, "fits"), 4L)
  expect_length(library, 12L)
  grid <- attr(library, "grid")
  expect_equal(nrow(grid), 12L)
  expect_identical(unique(grid$q), c(1L, 2L))
  expect_identical(unique(grid$kappa), c(0.05, 2))
  expect_identical(unique(grid$h), c(1L, 2L, 4L))
  expect_true(isTRUE(attr(library, "frozen")))
  expect_identical(
    attr(library, "tie_rule"),
    "first exact minimum in frozen candidate order"
  )
  expect_true(all(vapply(
    library,
    function(x) identical(x$source_periods, fixture$plan$fitting$periods),
    logical(1)
  )))
})

test_that("candidate construction is insulated from validation outcomes", {
  fixture <- make_workflow_fixture()
  altered_panel <- fixture$panel

  validation_dates <- fixture$plan$validation$raw_periods
  is_validation <- altered_panel$period %in% validation_dates
  altered_panel$outcome[is_validation] <-
    altered_panel$outcome[is_validation] + 1000 * seq_len(sum(is_validation))

  altered <- trsc_data(
    altered_panel,
    "unit", "period", "outcome",
    "Target", c("A", "B"),
    intervention_time = 16,
    pre_period = 1:15,
    predictors = "outcome"
  )
  altered_plan <- trsc_sample_plan(
    altered,
    calibration_periods = 1:4,
    fitting_periods = 6:9,
    validation_periods = 11:14,
    calibration_thin = 2,
    validation_thin = 2,
    target_period = 15,
    min_gap0 = 2,
    min_gap1 = 2
  )

  lib1 <- trsc_candidate_library(
    fixture$prepared, fixture$plan,
    q = c(1, 2), kappa = c(0.1, 1), bandwidths = c(1, 4)
  )
  lib2 <- trsc_candidate_library(
    altered, altered_plan,
    q = c(1, 2), kappa = c(0.1, 1), bandwidths = c(1, 4)
  )

  expect_equal(attr(lib1, "grid"), attr(lib2, "grid"))
  expect_equal(attr(lib1, "geometry")$Gamma,
               attr(lib2, "geometry")$Gamma, tolerance = 0)
  for (i in seq_along(lib1)) {
    expect_equal(
      as.numeric(lib1[[i]]$representative),
      as.numeric(lib2[[i]]$representative),
      tolerance = 1e-10
    )
    expect_equal(
      as.numeric(lib1[[i]]$metric_coordinates),
      as.numeric(lib2[[i]]$metric_coordinates),
      tolerance = 1e-10
    )
  }
})

test_that("calibration target outcomes are not used by geometry or fitting", {
  fixture <- make_workflow_fixture()
  altered_panel <- fixture$panel

  is_calibration_target <-
    altered_panel$unit == "Target" &
    altered_panel$period %in% fixture$plan$calibration$raw_periods
  altered_panel$outcome[is_calibration_target] <-
    altered_panel$outcome[is_calibration_target] + 5000

  altered <- trsc_data(
    altered_panel,
    "unit", "period", "outcome",
    "Target", c("A", "B"),
    intervention_time = 16,
    pre_period = 1:15,
    predictors = "outcome"
  )
  altered_plan <- trsc_sample_plan(
    altered, 1:4, 6:9, 11:14,
    calibration_thin = 2,
    validation_thin = 2,
    target_period = 15,
    min_gap0 = 2,
    min_gap1 = 2
  )

  lib1 <- trsc_candidate_library(
    fixture$prepared, fixture$plan,
    q = 1, kappa = 0.1, bandwidths = c(1, 4)
  )
  lib2 <- trsc_candidate_library(
    altered, altered_plan,
    q = 1, kappa = 0.1, bandwidths = c(1, 4)
  )

  expect_equal(attr(lib1, "geometry")$Gamma,
               attr(lib2, "geometry")$Gamma, tolerance = 0)
  expect_equal(
    as.numeric(lib1[[1]]$representative),
    as.numeric(lib2[[1]]$representative),
    tolerance = 1e-10
  )
})

test_that("end-to-end workflow freezes before validation and never refits", {
  fixture <- make_workflow_fixture()
  beta_assumption <- function(k) exp(-k / 2)

  workflow <- trsc_workflow(
    fixture$prepared,
    fixture$plan,
    q = c(1, 2),
    kappa = c(0.05, 1),
    bandwidths = c(1, 2, 4),
    mixing_beta = beta_assumption
  )

  expect_s3_class(workflow, "trsc_workflow")
  expect_s3_class(workflow$selection, "trsc_selection")
  expect_s3_class(workflow$selected, "trsc_readout")
  expect_s3_class(workflow$mixing, "trsc_mixing_bound")
  expect_identical(
    workflow$selection$validation_periods,
    fixture$plan$validation$retained_periods
  )
  expect_true(workflow$specification$candidate_order_fixed_before_validation)
  expect_false(workflow$specification$rolling_refit)
  expect_false(workflow$specification$final_full_sample_refit)
  expect_false(workflow$diagnostics$calibration_target_outcomes_used)
  expect_false(workflow$diagnostics$candidate_validation_outcomes_used)
  expect_true(workflow$diagnostics$validation_used_after_freeze)
  expect_false(workflow$diagnostics$theorem_assumptions_certified)
  expect_false(workflow$diagnostics$causal_identification)
  expect_identical(
    workflow$selected$source_periods,
    fixture$plan$fitting$periods
  )
})

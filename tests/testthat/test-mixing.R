make_mixing_plan <- function() {
  periods <- 1:20
  panel <- expand.grid(
    unit = c("Target", "A", "B"),
    period = periods,
    stringsAsFactors = FALSE
  )
  panel$outcome <- seq_len(nrow(panel))
  prepared <- trsc_data(
    panel,
    "unit", "period", "outcome",
    "Target", c("A", "B"),
    intervention_time = 21,
    pre_period = periods,
    predictors = "outcome"
  )
  plan <- trsc_sample_plan(
    prepared,
    calibration_periods = 1:6,
    fitting_periods = 9:13,
    validation_periods = 16:19,
    calibration_thin = 2,
    validation_thin = 2,
    target_period = 20,
    min_gap0 = 2,
    min_gap1 = 2
  )
  list(prepared = prepared, plan = plan)
}

test_that("mixing bound reproduces the four-term coupling expression", {
  fixture <- make_mixing_plan()
  beta <- function(k) 0.5^k
  bound <- trsc_mixing_bound(fixture$plan, beta)

  NC <- fixture$plan$calibration$retained_count
  v <- fixture$plan$validation$retained_count
  bC <- fixture$plan$specification$calibration_thin
  bV <- fixture$plan$specification$validation_thin
  g0 <- fixture$plan$gaps$g0_steps
  g1 <- fixture$plan$gaps$g1_steps

  expected <- (NC - 1) * beta(bC) + beta(g0) + beta(g1) +
    (v - 1) * beta(bV)

  expect_s3_class(bound, "trsc_mixing_bound")
  expect_equal(bound$delta_mix, expected, tolerance = 1e-15)
  expect_equal(
    sum(bound$contributions),
    bound$delta_mix,
    tolerance = 1e-15
  )
  expect_false(bound$certification$beta_function_verified)
  expect_false(bound$certification$stationarity_verified)
  expect_false(bound$certification$absolute_regularity_verified)
})

test_that("mixing bound rejects invalid assumed beta coefficients", {
  fixture <- make_mixing_plan()

  expect_error(
    trsc_mixing_bound(fixture$plan, function(k) -0.1),
    "finite value in \[0, 1\]"
  )
  expect_error(
    trsc_mixing_bound(fixture$plan, function(k) 1.2),
    "finite value in \[0, 1\]"
  )
})

test_that("explicit exponential-mixing design implements the paper formula", {
  design <- trsc_geometric_mixing_design(
    calibration_raw_length = 500,
    validation_raw_length = 200,
    C_beta = 1.5,
    tau_beta = 3,
    delta_mix = 0.05
  )

  expected_bC <- max(
    1L,
    as.integer(ceiling(3 * log(4 * 1.5 * 500 / 0.05)))
  )
  expected_bV <- max(
    1L,
    as.integer(ceiling(3 * log(4 * 1.5 * 200 / 0.05)))
  )
  expected_gap <- max(
    1L,
    as.integer(ceiling(3 * log(4 * 1.5 / 0.05)))
  )

  expect_s3_class(design, "trsc_geometric_mixing_design")
  expect_identical(design$calibration_thin, expected_bC)
  expect_identical(design$validation_thin, expected_bV)
  expect_identical(design$min_gap0, expected_gap)
  expect_identical(design$min_gap1, expected_gap)
  expect_lte(design$envelope_upper_bound, design$delta_mix_target)
  expect_true(design$target_respected_numerically)
  expect_false(design$assumption$certified)
})

test_that("explicit exponential-mixing design rejects an unusable validation segment", {
  expect_error(
    trsc_geometric_mixing_design(
      calibration_raw_length = 100,
      validation_raw_length = 3,
      C_beta = 1,
      tau_beta = 5,
      delta_mix = 0.01
    ),
    "retains fewer than two dates"
  )
})

test_that("TRSC data preparation and split contracts remain chronological", {
  panel <- expand.grid(
    unit = c("Target", "DonorA", "DonorB"),
    period = 2000:2007,
    stringsAsFactors = FALSE
  )
  panel$outcome <- seq_len(nrow(panel))
  panel$predictor <- 100 + seq_len(nrow(panel))

  prepared <- trsc_data(
    panel, "unit", "period", "outcome", "Target",
    c("DonorB", "DonorA"), 2008, 2000:2007,
    c("outcome", "predictor")
  )

  expect_s3_class(prepared, "trsc_data")
  expect_identical(prepared$schema_version, 2L)
  expect_s3_class(prepared$blocks, "trsc_blocks")
  expect_identical(prepared$blocks$schema_version, 2L)

  splits <- trsc_split(prepared, initial = 4, horizon = 2)
  expect_s3_class(splits, "trsc_splits")
  expect_identical(splits$schema_version, 1L)
  expect_identical(splits$folds[[1]]$train_period, 2000:2003)
  expect_identical(splits$folds[[1]]$validation_period, 2004:2005)
})

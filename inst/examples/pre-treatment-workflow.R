# Executable synthetic example for the stable pre-treatment TRSC workflow.

library(TRSC)

periods <- 1:15
A <- rep(0, length(periods))
B <- rep(10, length(periods))
target <- 10 * c(
  rep(0.2, 4), 0.2,
  0.2, 0.2, 0.2, 0.8,
  0.8, rep(0.8, 5)
)

panel <- expand.grid(
  unit = c("Target", "A", "B"),
  period = periods,
  stringsAsFactors = FALSE
)
panel$outcome <- c(rbind(target, A, B))

prepared <- trsc_data(
  panel, "unit", "period", "outcome",
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

workflow <- trsc_workflow(
  prepared,
  plan,
  q = c(1, 2),
  kappa = c(0.05, 1),
  bandwidths = c(1, 2, 4)
)

stopifnot(
  inherits(workflow, "trsc_workflow"),
  inherits(workflow$selected, "trsc_readout"),
  isFALSE(workflow$specification$rolling_refit),
  isFALSE(workflow$specification$final_full_sample_refit)
)

print(workflow)

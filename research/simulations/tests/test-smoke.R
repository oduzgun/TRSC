library(TRSC)

source("research/simulations/R/scenarios.R")
source("research/simulations/R/metrics.R")
source("research/simulations/R/benchmarks.R")
source("research/simulations/R/run-study.R")

required <- c(
  "constant", "affine", "curved", "abrupt",
  "exact_duplication", "near_redundancy",
  "weak_excitation", "reference_mismatch"
)

for (name in required) {
  scenario <- trsc_sim_scenario(
    scenario = name,
    T = 60,
    J = 5,
    seed = 101,
    noise_sd = 0.1,
    reference_size = 300
  )
  stopifnot(
    inherits(scenario, "trsc_sim_scenario"),
    identical(scenario$scenario, name),
    nrow(scenario$panel) == scenario$T * (scenario$J + 1L),
    nrow(scenario$reference_states) == 300
  )
}

dup <- trsc_sim_run_one(
  scenario_name = "exact_duplication",
  seed = 11,
  q = c(1L, 2L),
  kappa = c(0.1, 1),
  bandwidths = c(1L, 2L),
  T = 60,
  J = 5,
  noise_sd = 0.1,
  reference_size = 400,
  raw_benchmark_kappa = 0.3
)

stopifnot(
  inherits(dup$workflow, "trsc_workflow"),
  is.finite(dup$selected$cutoff_profile_loss),
  is.na(dup$selected$cutoff_raw_share_mse),
  all(is.na(dup$fit_metrics$raw_share_mse)),
  all(is.na(dup$benchmark_metrics$raw_share_mse)),
  isFALSE(dup$metadata$post_treatment_used),
  isFALSE(dup$metadata$causal_effect_claim)
)

smooth <- trsc_sim_run_one(
  scenario_name = "curved",
  seed = 12,
  q = c(1L, 2L),
  kappa = c(0.1, 1),
  bandwidths = c(1L, 2L),
  T = 60,
  J = 5,
  noise_sd = 0.1,
  reference_size = 400,
  raw_benchmark_kappa = 0.3
)

stopifnot(
  is.finite(smooth$selected$cutoff_profile_loss),
  is.finite(smooth$selected$cutoff_raw_share_mse),
  all(is.finite(smooth$fit_metrics$integrated_profile_loss)),
  all(is.finite(smooth$benchmark_metrics$integrated_profile_loss))
)

cat("TRSC simulation harness smoke checks passed.\n")

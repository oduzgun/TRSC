library(TRSC)

source("research/simulations/R/scenarios.R")
source("research/simulations/R/metrics.R")
source("research/simulations/R/benchmarks.R")
source("research/simulations/R/run-study.R")
source("research/simulations/R/reporting.R")

required <- c(
  "constant", "affine", "curved", "abrupt",
  "exact_duplication", "near_redundancy",
  "weak_excitation", "reference_mismatch",
  "interior", "persistent_errors"
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
    nrow(scenario$reference_states) == 300,
    nrow(scenario$heldout_states) == 300,
    identical(dim(scenario$Gamma_reference), c(scenario$J, scenario$J)),
    identical(dim(scenario$Gamma_heldout), c(scenario$J, scenario$J))
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
  is.finite(dup$selected$cutoff_heldout_prediction_loss),
  is.na(dup$selected$cutoff_raw_share_mse),
  all(is.na(dup$fit_metrics$raw_share_mse)),
  all(is.na(dup$benchmark_metrics$raw_share_mse)),
  nrow(dup$benchmark_cutoff) == 6L,
  all(is.finite(dup$benchmark_cutoff$cutoff_profile_loss)),
  all(is.finite(dup$benchmark_cutoff$cutoff_heldout_prediction_loss)),
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
  all(is.finite(smooth$fit_metrics$integrated_heldout_prediction_loss)),
  all(is.finite(smooth$benchmark_metrics$integrated_profile_loss)),
  all(is.finite(smooth$benchmark_metrics$integrated_heldout_prediction_loss)),
  nrow(smooth$benchmark_cutoff) == 6L
)

mismatch <- trsc_sim_scenario(
  scenario = "reference_mismatch",
  T = 60,
  J = 5,
  seed = 33,
  noise_sd = 0.1,
  reference_size = 500
)
stopifnot(
  max(abs(mismatch$Gamma_reference - mismatch$Gamma_heldout)) > 0.1
)

failed_grid <- trsc_sim_run_grid(
  scenarios = "constant",
  seeds = 1:2,
  T = 50,
  J = 5,
  reference_size = 300
)
stopifnot(
  nrow(failed_grid$failures) == 2L,
  failed_grid$metadata$attempted_runs == 2L,
  failed_grid$metadata$successful_runs == 0L,
  failed_grid$metadata$failed_runs == 2L,
  nrow(failed_grid$selected) == 0L
)

cat("TRSC simulation harness smoke checks passed.\n")


# Structural boundary scenarios are checked separately because the theorem-face
# construction can change the realized donor count.
separated <- trsc_sim_scenario(
  scenario = "separated_face",
  T = 60,
  J = 6,
  seed = 41,
  noise_sd = 0.05,
  reference_size = 400
)
vanishing <- trsc_sim_scenario(
  scenario = "vanishing_margin",
  T = 60,
  J = 6,
  seed = 41,
  noise_sd = 0.05,
  reference_size = 400
)
face_change <- trsc_sim_scenario(
  scenario = "face_change",
  T = 60,
  J = 6,
  seed = 42,
  noise_sd = 0.05,
  reference_size = 400
)
interior <- trsc_sim_scenario(
  scenario = "interior",
  T = 60,
  J = 6,
  seed = 43,
  noise_sd = 0.05,
  reference_size = 400
)

stopifnot(
  separated$structural$positive_normal_margin,
  separated$structural$theorem_class,
  !vanishing$structural$positive_normal_margin,
  !vanishing$structural$theorem_class,
  vanishing$structural$face_margin < separated$structural$face_margin,
  all(colSums(separated$true_weights) == 1),
  all(separated$true_weights[-c(1, 2), ] == 0),
  face_change$structural$face_before != face_change$structural$face_after,
  all(face_change$true_weights[1:2, 1:face_change$structural$switch_period] >= 0),
  all(face_change$true_weights[3:4, 1:face_change$structural$switch_period] == 0),
  all(face_change$true_weights[1:2, (face_change$structural$switch_period + 1):60] == 0),
  min(interior$true_weights) > 0
)

structural_run <- trsc_sim_run_one(
  scenario_name = "separated_face",
  seed = 44,
  q = c(1L, 2L),
  kappa = c(0.1, 1),
  bandwidths = c(1L, 2L),
  T = 60,
  J = 6,
  noise_sd = 0.05,
  reference_size = 400,
  raw_benchmark_kappa = 0.3,
  local_bandwidth = 4,
  tv_lambda = 0.3
)

stopifnot(
  nrow(structural_run$benchmark_metrics) == 6L,
  nrow(structural_run$benchmark_cutoff) == 6L,
  all(c("local_varying_simplex", "total_variation") %in%
        structural_run$benchmark_metrics$method),
  nrow(structural_run$selection_comparison) == 2L,
  identical(
    structural_run$selection_comparison$method,
    c("cutoff_targeted_local_linear", "ordinary_historical_validation")
  ),
  nrow(structural_run$design_diagnostics) == 1L,
  is.finite(structural_run$design_diagnostics$fitting_gram_min_eigenvalue),
  is.finite(structural_run$design_diagnostics$reference_min_eigenvalue)
)

weak_run <- trsc_sim_run_one(
  scenario_name = "weak_excitation",
  seed = 45,
  q = 1L,
  kappa = c(0.1, 1),
  bandwidths = c(1L, 2L),
  T = 60,
  J = 6,
  noise_sd = 0.05,
  reference_size = 400,
  local_bandwidth = 4,
  tv_lambda = 0.3
)
stopifnot(
  weak_run$design_diagnostics$contrast_min_eigenvalue <
    weak_run$design_diagnostics$reference_min_eigenvalue
)

# Direct solver sanity checks for the added research benchmarks.
objects <- trsc_sim_prepare(interior)
local_fit <- trsc_sim_local_varying_simplex(
  objects$prepared, objects$plan$fitting$periods, bandwidth = 4
)
tv_fit <- trsc_sim_total_variation_fit(
  objects$prepared, objects$plan$fitting$periods, lambda = 0.3
)
stopifnot(
  max(abs(colSums(local_fit$weights) - 1)) < 1e-6,
  min(local_fit$weights) > -1e-6,
  max(abs(colSums(tv_fit$weights) - 1)) < 1e-6,
  min(tv_fit$weights) > -1e-6
)


report_grid <- trsc_sim_run_grid(
  scenarios = "constant",
  seeds = 71:72,
  q = 1L,
  kappa = 0.2,
  bandwidths = 1L,
  T = 60,
  J = 4,
  noise_sd = 0.05,
  reference_size = 250,
  raw_benchmark_kappa = 0.2,
  local_bandwidth = 4,
  tv_lambda = 0.2
)

stopifnot(
  report_grid$metadata$attempted_runs == 2L,
  report_grid$metadata$successful_runs == 2L,
  report_grid$metadata$failed_runs == 0L,
  all(c("scenario", "seed") %in% names(report_grid$selection_comparison)),
  all(c("scenario", "seed") %in% names(report_grid$design_diagnostics)),
  all(c("scenario", "seed") %in% names(report_grid$solver_diagnostics))
)

report <- trsc_sim_mc_report(report_grid)
stopifnot(
  inherits(report, "trsc_sim_mc_report"),
  nrow(report$failures) == 1L,
  report$failures$accounting_complete,
  report$failures$failure_rate == 0,
  nrow(report$selected_cutoff) > 0L,
  all(report$selected_cutoff$n == 2L |
        report$selected_cutoff$metric == "cutoff_raw_share_mse"),
  nrow(report$selection_frequency) > 0L,
  isFALSE(report$manifest$simulation_findings_claimed),
  !is.na(report$manifest$package_version)
)

frequency_sums <- aggregate(
  proportion ~ scenario + method,
  data = report$selection_frequency,
  FUN = sum
)
stopifnot(max(abs(frequency_sums$proportion - 1)) < 1e-12)

report_dir <- tempfile("trsc-report-")
trsc_sim_write_report_bundle(report, report_dir)
stopifnot(
  file.exists(file.path(report_dir, "selected_cutoff.csv")),
  file.exists(file.path(report_dir, "solver.csv")),
  file.exists(file.path(report_dir, "failures.csv")),
  file.exists(file.path(report_dir, "reproducibility_manifest.rds")),
  file.exists(file.path(report_dir, "sessionInfo.txt"))
)
unlink(report_dir, recursive = TRUE)

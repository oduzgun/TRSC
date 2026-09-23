# Research-only Monte Carlo orchestration.
#
# Source scenarios.R, metrics.R and benchmarks.R before using these functions.

trsc_sim_prepare <- function(scenario,
                             calibration_thin = 2L,
                             validation_thin = 2L) {
  if (!inherits(scenario, "trsc_sim_scenario")) {
    stop("scenario must be a trsc_sim_scenario object.", call. = FALSE)
  }

  prepared <- trsc_data(
    scenario$panel,
    unit = "unit",
    time = "period",
    outcome = "outcome",
    treated = "Target",
    donors = scenario$donor_names,
    intervention_time = scenario$T + 1L,
    pre_period = seq_len(scenario$T),
    predictors = "outcome"
  )

  index <- scenario$plan_index
  plan <- trsc_sample_plan(
    prepared,
    calibration_periods = index$calibration,
    fitting_periods = index$fitting,
    validation_periods = index$validation,
    calibration_thin = calibration_thin,
    validation_thin = validation_thin,
    target_period = index$target,
    min_gap0 = 2L,
    min_gap1 = 2L
  )
  list(prepared = prepared, plan = plan)
}

trsc_sim_run_one <- function(
    scenario_name,
    seed,
    q = c(1L, 2L),
    kappa = c(0.03, 0.3, 3),
    bandwidths = NULL,
    T = 90L,
    J = 6L,
    noise_sd = 0.20,
    reference_size = 4000L,
    raw_benchmark_kappa = 0.3) {

  scenario <- trsc_sim_scenario(
    scenario = scenario_name,
    T = T,
    J = J,
    seed = seed,
    noise_sd = noise_sd,
    reference_size = reference_size
  )
  objects <- trsc_sim_prepare(scenario)

  workflow <- trsc_workflow(
    objects$prepared,
    objects$plan,
    q = q,
    kappa = kappa,
    bandwidths = bandwidths
  )
  workflow_metrics <- trsc_sim_evaluate_workflow(workflow, scenario)
  benchmark <- trsc_sim_benchmark_metrics(
    objects$prepared,
    scenario,
    objects$plan,
    raw_kappa = raw_benchmark_kappa
  )

  selected_row <- data.frame(
    scenario = scenario$scenario,
    seed = scenario$seed,
    path_type = scenario$path_type,
    representation_identified = scenario$representation_identified,
    method = "TRSC_selected",
    selected_q = workflow_metrics$cutoff$selected_q,
    selected_kappa = workflow_metrics$cutoff$selected_kappa,
    selected_h = workflow_metrics$cutoff$selected_h,
    cutoff_profile_loss = workflow_metrics$cutoff$cutoff_profile_loss,
    cutoff_heldout_prediction_loss =
      workflow_metrics$cutoff$cutoff_heldout_prediction_loss,
    cutoff_raw_share_mse = workflow_metrics$cutoff$cutoff_raw_share_mse,
    validation_score = workflow_metrics$cutoff$validation_score,
    stringsAsFactors = FALSE
  )

  fit_rows <- workflow_metrics$fits
  fit_rows$scenario <- scenario$scenario
  fit_rows$seed <- scenario$seed
  fit_rows$path_type <- scenario$path_type
  fit_rows$representation_identified <- scenario$representation_identified
  fit_rows$method <- "TRSC_candidate"
  fit_rows <- fit_rows[, c(
    "scenario", "seed", "path_type", "representation_identified",
    "method", "fit", "q", "kappa",
    "integrated_profile_loss", "integrated_heldout_prediction_loss",
    "raw_share_mse", "training_rmse", "temporal_roughness"
  )]

  benchmark_path <- benchmark$path
  benchmark_path$scenario <- scenario$scenario
  benchmark_path$seed <- scenario$seed
  benchmark_path$path_type <- scenario$path_type
  benchmark_path$representation_identified <-
    scenario$representation_identified
  benchmark_path <- benchmark_path[, c(
    "scenario", "seed", "path_type", "representation_identified",
    "method", "integrated_profile_loss", "integrated_heldout_prediction_loss",
    "raw_share_mse", "training_rmse"
  )]

  benchmark_cutoff <- benchmark$cutoff
  benchmark_cutoff$scenario <- scenario$scenario
  benchmark_cutoff$seed <- scenario$seed
  benchmark_cutoff$path_type <- scenario$path_type
  benchmark_cutoff$representation_identified <-
    scenario$representation_identified

  list(
    scenario = scenario,
    prepared = objects$prepared,
    plan = objects$plan,
    workflow = workflow,
    selected = selected_row,
    fit_metrics = fit_rows,
    benchmark_metrics = benchmark_path,
    benchmark_cutoff = benchmark_cutoff,
    metadata = list(
      generated_at = "not recorded: deterministic research object",
      q = q,
      kappa = kappa,
      bandwidths = if (is.null(bandwidths)) "default_dyadic_plus_full"
                   else bandwidths,
      raw_benchmark_kappa = raw_benchmark_kappa,
      post_treatment_used = FALSE,
      causal_effect_claim = FALSE
    )
  )
}

trsc_sim_run_grid <- function(
    scenarios = c(
      "constant", "affine", "curved", "abrupt",
      "exact_duplication", "near_redundancy",
      "weak_excitation", "reference_mismatch"
    ),
    seeds = 1:20,
    output_dir = NULL,
    ...) {

  if (!is.character(scenarios) || !length(scenarios) ||
      anyNA(scenarios) || any(!nzchar(scenarios))) {
    stop("scenarios must be a nonempty character vector.", call. = FALSE)
  }
  if (!is.numeric(seeds) || !length(seeds) || anyNA(seeds) ||
      any(!is.finite(seeds))) {
    stop("seeds must be a nonempty finite numeric vector.", call. = FALSE)
  }

  runs <- vector("list", length(scenarios) * length(seeds))
  failures <- list()
  index <- 0L
  failure_index <- 0L
  for (scenario_name in scenarios) {
    for (seed in seeds) {
      index <- index + 1L
      attempt <- tryCatch(
        trsc_sim_run_one(
          scenario_name = scenario_name,
          seed = as.integer(seed),
          ...
        ),
        error = function(e) e
      )
      if (inherits(attempt, "error")) {
        failure_index <- failure_index + 1L
        failures[[failure_index]] <- data.frame(
          scenario = scenario_name,
          seed = as.integer(seed),
          error_class = class(attempt)[1L],
          message = conditionMessage(attempt),
          stringsAsFactors = FALSE
        )
        runs[index] <- list(NULL)
      } else {
        runs[index] <- list(attempt)
      }
    }
  }

  successful_runs <- Filter(Negate(is.null), runs)
  bind_or_empty <- function(component) {
    if (!length(successful_runs)) return(data.frame())
    out <- do.call(rbind, lapply(successful_runs, function(x) x[[component]]))
    rownames(out) <- NULL
    out
  }

  selected <- bind_or_empty("selected")
  fit_metrics <- bind_or_empty("fit_metrics")
  benchmark_metrics <- bind_or_empty("benchmark_metrics")
  benchmark_cutoff <- bind_or_empty("benchmark_cutoff")
  failure_table <- if (length(failures)) {
    out <- do.call(rbind, failures)
    rownames(out) <- NULL
    out
  } else {
    data.frame(
      scenario = character(),
      seed = integer(),
      error_class = character(),
      message = character(),
      stringsAsFactors = FALSE
    )
  }

  result <- list(
    selected = selected,
    fit_metrics = fit_metrics,
    benchmark_metrics = benchmark_metrics,
    benchmark_cutoff = benchmark_cutoff,
    failures = failure_table,
    runs = runs,
    metadata = list(
      scenarios = scenarios,
      seeds = as.integer(seeds),
      attempted_runs = length(scenarios) * length(seeds),
      successful_runs = length(successful_runs),
      failed_runs = nrow(failure_table),
      simulation_findings_claimed = FALSE,
      purpose = "research design and software discrimination"
    )
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(
      selected,
      file.path(output_dir, "selected_cutoff_metrics.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      fit_metrics,
      file.path(output_dir, "candidate_path_metrics.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      benchmark_metrics,
      file.path(output_dir, "benchmark_path_metrics.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      benchmark_cutoff,
      file.path(output_dir, "benchmark_cutoff_metrics.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      failure_table,
      file.path(output_dir, "simulation_failures.csv"),
      row.names = FALSE
    )
    saveRDS(result$metadata, file.path(output_dir, "simulation_metadata.rds"))
  }

  result
}

# Research-only reporting helpers for TRSC simulation studies.
#
# These functions summarize completed Monte Carlo objects. They do not generate
# inferential claims automatically and are intentionally outside the package
# namespace.

.trsc_sim_validate_grid_result <- function(result) {
  required <- c(
    "selected", "fit_metrics", "benchmark_metrics", "benchmark_cutoff",
    "selection_comparison", "design_diagnostics", "solver_diagnostics",
    "failures", "metadata"
  )
  if (!is.list(result) || !all(required %in% names(result))) {
    stop("result must be a trsc_sim_run_grid() result.", call. = FALSE)
  }
  invisible(result)
}

.trsc_sim_numeric_summary <- function(data, groups, metrics,
                                      conf_level = 0.95) {
  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  if (!length(metrics)) return(data.frame())
  if (!all(groups %in% names(data)) || !all(metrics %in% names(data))) {
    stop("Requested summary columns are unavailable.", call. = FALSE)
  }
  if (length(conf_level) != 1L || !is.numeric(conf_level) ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (!nrow(data)) return(data.frame())

  key <- interaction(data[groups], drop = TRUE, lex.order = TRUE)
  index_groups <- split(seq_len(nrow(data)), key, drop = TRUE)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  rows <- list()
  counter <- 0L
  for (idx in index_groups) {
    base <- data[idx[1L], groups, drop = FALSE]
    for (metric in metrics) {
      x <- data[[metric]][idx]
      usable <- is.finite(x)
      values <- x[usable]
      counter <- counter + 1L
      if (!length(values)) {
        stats_row <- data.frame(
          metric = metric,
          n = 0L,
          n_missing = length(x),
          mean = NA_real_,
          sd = NA_real_,
          mcse = NA_real_,
          median = NA_real_,
          empirical_q025 = NA_real_,
          empirical_q975 = NA_real_,
          mean_ci_lower = NA_real_,
          mean_ci_upper = NA_real_,
          stringsAsFactors = FALSE
        )
      } else {
        n <- length(values)
        sdev <- if (n > 1L) stats::sd(values) else NA_real_
        mcse <- if (n > 1L) sdev / sqrt(n) else NA_real_
        mean_value <- mean(values)
        quant <- stats::quantile(
          values, probs = c(0.025, 0.975),
          names = FALSE, type = 8
        )
        stats_row <- data.frame(
          metric = metric,
          n = n,
          n_missing = length(x) - n,
          mean = mean_value,
          sd = sdev,
          mcse = mcse,
          median = stats::median(values),
          empirical_q025 = quant[1L],
          empirical_q975 = quant[2L],
          mean_ci_lower = if (is.finite(mcse)) mean_value - z * mcse else NA_real_,
          mean_ci_upper = if (is.finite(mcse)) mean_value + z * mcse else NA_real_,
          stringsAsFactors = FALSE
        )
      }
      rows[[counter]] <- cbind(base, stats_row, stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

trsc_sim_failure_summary <- function(result) {
  .trsc_sim_validate_grid_result(result)
  scenarios <- as.character(result$metadata$scenarios)
  seeds <- as.integer(result$metadata$seeds)
  failures <- result$failures
  successes <- result$selected

  do.call(rbind, lapply(scenarios, function(scenario) {
    failed <- if (nrow(failures)) sum(failures$scenario == scenario) else 0L
    successful <- if (nrow(successes)) {
      sum(successes$scenario == scenario)
    } else {
      0L
    }
    attempted <- length(seeds)
    data.frame(
      scenario = scenario,
      attempted = attempted,
      successful = successful,
      failed = failed,
      failure_rate = if (attempted) failed / attempted else NA_real_,
      accounting_complete = identical(
        as.integer(successful + failed), as.integer(attempted)
      ),
      stringsAsFactors = FALSE
    )
  }))
}

trsc_sim_selection_frequency <- function(result) {
  .trsc_sim_validate_grid_result(result)
  x <- result$selection_comparison
  if (!nrow(x)) return(data.frame())
  cols <- c(
    "scenario", "method", "selected_q", "selected_kappa", "selected_h"
  )
  key <- interaction(x[cols], drop = TRUE, lex.order = TRUE)
  pieces <- split(seq_len(nrow(x)), key, drop = TRUE)
  out <- do.call(rbind, lapply(pieces, function(idx) {
    base <- x[idx[1L], cols, drop = FALSE]
    denominator <- sum(
      x$scenario == base$scenario[1L] &
        x$method == base$method[1L]
    )
    data.frame(
      base,
      count = length(idx),
      proportion = length(idx) / denominator,
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out[order(out$scenario, out$method, -out$proportion), , drop = FALSE]
}

trsc_sim_reproducibility_manifest <- function(result,
                                              git_sha = Sys.getenv(
                                                "GITHUB_SHA",
                                                unset = NA_character_
                                              )) {
  .trsc_sim_validate_grid_result(result)
  package_version <- tryCatch(
    as.character(utils::packageVersion("TRSC")),
    error = function(e) NA_character_
  )
  list(
    schema_version = 1L,
    package = "TRSC",
    package_version = package_version,
    git_sha = if (nzchar(git_sha)) git_sha else NA_character_,
    r_version = R.version.string,
    platform = R.version$platform,
    scenarios = as.character(result$metadata$scenarios),
    seeds = as.integer(result$metadata$seeds),
    attempted_runs = result$metadata$attempted_runs,
    successful_runs = result$metadata$successful_runs,
    failed_runs = result$metadata$failed_runs,
    simulation_findings_claimed = FALSE,
    session_info = capture.output(utils::sessionInfo())
  )
}

trsc_sim_mc_report <- function(result, conf_level = 0.95) {
  .trsc_sim_validate_grid_result(result)

  selected_summary <- .trsc_sim_numeric_summary(
    result$selected,
    groups = c("scenario", "method"),
    metrics = c(
      "cutoff_profile_loss",
      "cutoff_heldout_prediction_loss",
      "cutoff_raw_share_mse",
      "validation_score"
    ),
    conf_level = conf_level
  )

  selection_target_summary <- .trsc_sim_numeric_summary(
    result$selection_comparison,
    groups = c("scenario", "method"),
    metrics = c(
      "cutoff_profile_loss",
      "cutoff_heldout_prediction_loss",
      "cutoff_raw_share_mse",
      "validation_score"
    ),
    conf_level = conf_level
  )

  benchmark_path_summary <- .trsc_sim_numeric_summary(
    result$benchmark_metrics,
    groups = c("scenario", "method"),
    metrics = c(
      "integrated_profile_loss",
      "integrated_heldout_prediction_loss",
      "raw_share_mse",
      "training_rmse"
    ),
    conf_level = conf_level
  )

  benchmark_cutoff_summary <- .trsc_sim_numeric_summary(
    result$benchmark_cutoff,
    groups = c("scenario", "method"),
    metrics = c(
      "cutoff_profile_loss",
      "cutoff_heldout_prediction_loss",
      "cutoff_raw_share_mse"
    ),
    conf_level = conf_level
  )

  design_summary <- .trsc_sim_numeric_summary(
    result$design_diagnostics,
    groups = c("scenario"),
    metrics = c(
      "fitting_gram_rank",
      "fitting_gram_min_eigenvalue",
      "contrast_rank",
      "contrast_min_eigenvalue",
      "calibration_rank",
      "reference_rank",
      "reference_min_eigenvalue"
    ),
    conf_level = conf_level
  )

  solver_summary <- .trsc_sim_numeric_summary(
    result$solver_diagnostics,
    groups = c("scenario", "q", "kappa"),
    metrics = c(
      "iterations",
      "primal_residual",
      "dual_residual",
      "maximum_simplex_violation"
    ),
    conf_level = conf_level
  )

  structure(
    list(
      schema_version = 1L,
      conf_level = conf_level,
      selected_cutoff = selected_summary,
      selection_target_comparison = selection_target_summary,
      benchmark_path = benchmark_path_summary,
      benchmark_cutoff = benchmark_cutoff_summary,
      design = design_summary,
      solver = solver_summary,
      failures = trsc_sim_failure_summary(result),
      selection_frequency = trsc_sim_selection_frequency(result),
      manifest = trsc_sim_reproducibility_manifest(result),
      interpretation = list(
        monte_carlo_ci = "normal approximation for the simulation mean using the Monte Carlo standard error",
        empirical_interval = "2.5% and 97.5% empirical quantiles across successful replications",
        failures_retained = TRUE,
        causal_claim = FALSE,
        theorem_certification = FALSE
      )
    ),
    class = "trsc_sim_mc_report"
  )
}

trsc_sim_write_report_bundle <- function(report, output_dir) {
  if (!inherits(report, "trsc_sim_mc_report") || !is.list(report)) {
    stop("report must be a trsc_sim_mc_report object.", call. = FALSE)
  }
  if (length(output_dir) != 1L || !is.character(output_dir) ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be one nonempty path.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  tables <- c(
    "selected_cutoff", "selection_target_comparison",
    "benchmark_path", "benchmark_cutoff",
    "design", "solver", "failures", "selection_frequency"
  )
  for (name in tables) {
    utils::write.csv(
      report[[name]],
      file.path(output_dir, paste0(name, ".csv")),
      row.names = FALSE
    )
  }
  saveRDS(report$manifest, file.path(output_dir, "reproducibility_manifest.rds"))
  writeLines(
    report$manifest$session_info,
    con = file.path(output_dir, "sessionInfo.txt")
  )
  invisible(normalizePath(output_dir, mustWork = TRUE))
}

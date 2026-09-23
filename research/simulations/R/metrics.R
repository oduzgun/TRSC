# Research-only losses for TRSC simulations.

trsc_sim_profile_loss <- function(weights_estimated, weights_true, Gamma) {
  if (!is.matrix(Gamma) || nrow(Gamma) != ncol(Gamma)) {
    stop("Gamma must be square.", call. = FALSE)
  }
  if (is.vector(weights_estimated)) weights_estimated <- matrix(weights_estimated, ncol = 1L)
  if (is.vector(weights_true)) weights_true <- matrix(weights_true, ncol = 1L)
  if (!is.matrix(weights_estimated) || !is.matrix(weights_true) ||
      !identical(dim(weights_estimated), dim(weights_true)) ||
      nrow(weights_estimated) != nrow(Gamma)) {
    stop("Estimated and true weights must be aligned donor-by-date matrices.", call. = FALSE)
  }
  E <- t(weights_estimated - weights_true)
  mean(rowSums((E %*% Gamma) * E))
}

trsc_sim_reference_prediction_loss <- function(w_estimated, w_true, reference_states) {
  if (!is.numeric(w_estimated) || !is.numeric(w_true) ||
      length(w_estimated) != length(w_true) ||
      !is.matrix(reference_states) ||
      ncol(reference_states) != length(w_true)) {
    stop("Weights and reference states are not aligned.", call. = FALSE)
  }
  diff <- as.numeric(w_estimated - w_true)
  mean(drop(reference_states %*% diff)^2)
}

trsc_sim_share_loss <- function(weights_estimated, weights_true,
                                identified = FALSE) {
  if (!isTRUE(identified)) return(NA_real_)
  if (is.vector(weights_estimated)) weights_estimated <- matrix(weights_estimated, ncol = 1L)
  if (is.vector(weights_true)) weights_true <- matrix(weights_true, ncol = 1L)
  if (!identical(dim(weights_estimated), dim(weights_true))) {
    stop("Estimated and true weights must have the same dimensions.", call. = FALSE)
  }
  mean((weights_estimated - weights_true)^2)
}

trsc_sim_evaluate_workflow <- function(workflow, scenario) {
  if (!inherits(workflow, "trsc_workflow")) {
    stop("workflow must be a trsc_workflow object.", call. = FALSE)
  }
  if (!inherits(scenario, "trsc_sim_scenario")) {
    stop("scenario must be a trsc_sim_scenario object.", call. = FALSE)
  }

  fit_periods <- workflow$plan$fitting$periods
  true_fit <- scenario$true_weights[, fit_periods, drop = FALSE]
  fits <- attr(workflow$candidate_library, "fits")

  fit_metrics <- do.call(rbind, lapply(names(fits), function(name) {
    fit <- fits[[name]]
    estimated <- unclass(stats::weights(fit))
    data.frame(
      fit = name,
      q = fit$specification$q,
      kappa = fit$specification$kappa,
      integrated_profile_loss = trsc_sim_profile_loss(
        estimated, true_fit, scenario$Gamma_reference
      ),
      integrated_heldout_prediction_loss = trsc_sim_profile_loss(
        estimated, true_fit, scenario$Gamma_heldout
      ),
      raw_share_mse = trsc_sim_share_loss(
        estimated, true_fit, scenario$representation_identified
      ),
      training_rmse = sqrt(mean(stats::residuals(fit)^2)),
      temporal_roughness = fit$objective$temporal_roughness,
      stringsAsFactors = FALSE
    )
  }))

  target_period <- workflow$plan$target_period
  true_cutoff <- scenario$true_weights[, target_period]
  selected <- workflow$selected$representative[scenario$donor_names]

  cutoff <- data.frame(
    selected_candidate = workflow$selection$selected_name,
    selected_q = workflow$selected$candidate_specification$q,
    selected_kappa = workflow$selected$candidate_specification$kappa,
    selected_h = workflow$selected$candidate_specification$h,
    cutoff_profile_loss = trsc_sim_reference_prediction_loss(
      selected, true_cutoff, scenario$reference_states
    ),
    cutoff_heldout_prediction_loss = trsc_sim_reference_prediction_loss(
      selected, true_cutoff, scenario$heldout_states
    ),
    cutoff_raw_share_mse = if (scenario$representation_identified) {
      mean((selected - true_cutoff)^2)
    } else {
      NA_real_
    },
    validation_score = workflow$selection$scores[
      workflow$selection$selected_name
    ],
    stringsAsFactors = FALSE
  )

  list(fits = fit_metrics, cutoff = cutoff)
}


trsc_sim_selection_comparison <- function(prepared, workflow, scenario) {
  if (!inherits(prepared, "trsc_data") ||
      !inherits(workflow, "trsc_workflow") ||
      !inherits(scenario, "trsc_sim_scenario")) {
    stop("prepared, workflow and scenario are not compatible.", call. = FALSE)
  }
  target_period <- workflow$plan$target_period
  true_cutoff <- scenario$true_weights[, target_period]

  historical <- trsc_sim_historical_select(
    prepared,
    workflow$candidate_library,
    workflow$plan$validation$retained_periods
  )

  evaluate <- function(label, candidate_name, candidate, score, target) {
    w <- candidate$representative[scenario$donor_names]
    data.frame(
      method = label,
      selected_candidate = candidate_name,
      selected_q = candidate$candidate_specification$q,
      selected_kappa = candidate$candidate_specification$kappa,
      selected_h = candidate$candidate_specification$h,
      validation_target = target,
      validation_score = as.numeric(score),
      cutoff_profile_loss = trsc_sim_reference_prediction_loss(
        w, true_cutoff, scenario$reference_states
      ),
      cutoff_heldout_prediction_loss = trsc_sim_reference_prediction_loss(
        w, true_cutoff, scenario$heldout_states
      ),
      cutoff_raw_share_mse = if (scenario$representation_identified) {
        mean((w - true_cutoff)^2)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }

  rbind(
    evaluate(
      "cutoff_targeted_local_linear",
      workflow$selection$selected_name,
      workflow$selected,
      workflow$selection$scores[workflow$selection$selected_name],
      "cutoff risk"
    ),
    evaluate(
      "ordinary_historical_validation",
      historical$selected_name,
      historical$selected,
      historical$scores[historical$selected_name],
      historical$target
    )
  )
}


trsc_sim_solver_diagnostics <- function(workflow, scenario) {
  if (!inherits(workflow, "trsc_workflow") ||
      !inherits(scenario, "trsc_sim_scenario")) {
    stop("workflow and scenario are not compatible.", call. = FALSE)
  }
  fits <- attr(workflow$candidate_library, "fits")
  out <- do.call(rbind, lapply(names(fits), function(name) {
    fit <- fits[[name]]
    data.frame(
      scenario = scenario$scenario,
      seed = scenario$seed,
      fit = name,
      q = fit$specification$q,
      kappa = fit$specification$kappa,
      solver_status = as.character(fit$solver$status),
      iterations = as.integer(fit$solver$iterations),
      primal_residual = as.numeric(fit$solver$primal_residual),
      dual_residual = as.numeric(fit$solver$dual_residual),
      maximum_simplex_violation =
        as.numeric(fit$diagnostics$maximum_simplex_violation),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

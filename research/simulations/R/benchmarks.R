# Research-only matched-information benchmarks for TRSC simulations.

trsc_sim_static_simplex <- function(prepared, fitting_periods,
                                    solver_control = list()) {
  if (!inherits(prepared, "trsc_data")) {
    stop("prepared must be a trsc_data object.", call. = FALSE)
  }
  index <- match(fitting_periods, prepared$blocks$periods)
  if (anyNA(index)) stop("fitting_periods are not available.", call. = FALSE)

  X <- prepared$blocks$Y0[index, , drop = FALSE]
  y <- as.numeric(prepared$blocks$Y1[index])
  J <- ncol(X)

  P <- 2 * Matrix::Matrix(crossprod(X), sparse = TRUE)
  P <- Matrix::forceSymmetric(P, uplo = "U")
  qvec <- -2 * as.numeric(crossprod(X, y))
  A <- rbind(
    Matrix::Matrix(matrix(1, nrow = 1L, ncol = J), sparse = TRUE),
    Matrix::Diagonal(J)
  )
  l <- c(1, rep(0, J))
  u <- c(1, rep(1, J))

  defaults <- list(
    verbose = FALSE,
    eps_abs = 1e-8,
    eps_rel = 1e-8,
    max_iter = 100000L,
    polishing = TRUE
  )
  settings <- do.call(
    osqp::osqpSettings,
    utils::modifyList(defaults, solver_control)
  )
  result <- osqp::solve_osqp(P = P, q = qvec, A = A, l = l, u = u,
                             pars = settings)
  if (tolower(trimws(result$info$status)) != "solved") {
    stop("Static benchmark failed with OSQP status: ", result$info$status,
         call. = FALSE)
  }

  weights <- as.numeric(result$x)
  names(weights) <- colnames(X)
  fitted <- drop(X %*% weights)

  structure(
    list(
      weights = weights,
      fitted = fitted,
      residuals = y - fitted,
      periods = fitting_periods,
      solver_status = result$info$status,
      information = "same fitting outcomes as TRSC"
    ),
    class = "trsc_sim_static_fit"
  )
}

trsc_sim_raw_share_fit <- function(prepared, fitting_periods, kappa, q) {
  donor_names <- colnames(prepared$blocks$Y0)
  J <- length(donor_names)
  identity_states <- sqrt(J) * diag(J)
  colnames(identity_states) <- donor_names
  raw_geometry <- trsc_geometry(identity_states)
  fit <- trsc_fit(
    prepared,
    geometry = raw_geometry,
    kappa = kappa,
    q = q,
    periods = fitting_periods
  )
  fit$research_benchmark <- list(
    penalty = "Euclidean raw-share differences",
    representation_invariant = FALSE
  )
  fit
}

trsc_sim_benchmark_metrics <- function(prepared, scenario, plan,
                                       raw_kappa = 1) {
  true_fit <- scenario$true_weights[, plan$fitting$periods, drop = FALSE]
  true_cutoff <- scenario$true_weights[, plan$target_period]

  static <- trsc_sim_static_simplex(prepared, plan$fitting$periods)
  static_path <- matrix(
    static$weights,
    nrow = length(static$weights),
    ncol = length(plan$fitting$periods),
    dimnames = list(names(static$weights), as.character(plan$fitting$periods))
  )
  static_cutoff <- static$weights

  calibration_geometry <- trsc_calibration_geometry(prepared, plan)
  period_specific <- trsc_fit(
    prepared,
    calibration_geometry,
    kappa = 0,
    q = 1,
    periods = plan$fitting$periods
  )

  raw_q1 <- trsc_sim_raw_share_fit(
    prepared, plan$fitting$periods, raw_kappa, 1
  )
  raw_q2 <- trsc_sim_raw_share_fit(
    prepared, plan$fitting$periods, raw_kappa, 2
  )

  summarize_path <- function(label, W, training_rmse) {
    data.frame(
      method = label,
      integrated_profile_loss = trsc_sim_profile_loss(
        W, true_fit, scenario$Gamma_reference
      ),
      raw_share_mse = trsc_sim_share_loss(
        W, true_fit, scenario$representation_identified
      ),
      training_rmse = training_rmse,
      stringsAsFactors = FALSE
    )
  }

  path_metrics <- rbind(
    summarize_path(
      "static_simplex",
      static_path,
      sqrt(mean(static$residuals^2))
    ),
    summarize_path(
      "period_specific",
      unclass(stats::weights(period_specific)),
      sqrt(mean(stats::residuals(period_specific)^2))
    ),
    summarize_path(
      "raw_share_q1",
      unclass(stats::weights(raw_q1)),
      sqrt(mean(stats::residuals(raw_q1)^2))
    ),
    summarize_path(
      "raw_share_q2",
      unclass(stats::weights(raw_q2)),
      sqrt(mean(stats::residuals(raw_q2)^2))
    )
  )

  cutoff_metrics <- data.frame(
    method = "static_simplex",
    cutoff_profile_loss = trsc_sim_reference_prediction_loss(
      static_cutoff, true_cutoff, scenario$reference_states
    ),
    cutoff_raw_share_mse = if (scenario$representation_identified) {
      mean((static_cutoff - true_cutoff)^2)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )

  list(
    path = path_metrics,
    cutoff = cutoff_metrics,
    objects = list(
      static = static,
      period_specific = period_specific,
      raw_q1 = raw_q1,
      raw_q2 = raw_q2
    )
  )
}

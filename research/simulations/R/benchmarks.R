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
                                       raw_kappa = 1,
                                       local_bandwidth = 5,
                                       tv_lambda = 1) {
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
  local_varying <- trsc_sim_local_varying_simplex(
    prepared, plan$fitting$periods, bandwidth = local_bandwidth
  )
  total_variation <- trsc_sim_total_variation_fit(
    prepared, plan$fitting$periods, lambda = tv_lambda
  )

  summarize_path <- function(label, W, training_rmse) {
    data.frame(
      method = label,
      integrated_profile_loss = trsc_sim_profile_loss(
        W, true_fit, scenario$Gamma_reference
      ),
      integrated_heldout_prediction_loss = trsc_sim_profile_loss(
        W, true_fit, scenario$Gamma_heldout
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
    ),
    summarize_path(
      "local_varying_simplex",
      local_varying$weights,
      sqrt(mean(local_varying$residuals^2))
    ),
    summarize_path(
      "total_variation",
      total_variation$weights,
      sqrt(mean(total_variation$residuals^2))
    )
  )

  endpoint_weights <- list(
    static_simplex = static_cutoff,
    period_specific = unclass(stats::weights(period_specific))[, ncol(period_specific$weights)],
    raw_share_q1 = unclass(stats::weights(raw_q1))[, ncol(raw_q1$weights)],
    raw_share_q2 = unclass(stats::weights(raw_q2))[, ncol(raw_q2$weights)],
    local_varying_simplex = local_varying$weights[, ncol(local_varying$weights)],
    total_variation = total_variation$weights[, ncol(total_variation$weights)]
  )

  cutoff_metrics <- do.call(rbind, lapply(names(endpoint_weights), function(label) {
    w <- endpoint_weights[[label]]
    data.frame(
      method = label,
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
  }))
  rownames(cutoff_metrics) <- NULL

  list(
    path = path_metrics,
    cutoff = cutoff_metrics,
    objects = list(
      static = static,
      period_specific = period_specific,
      raw_q1 = raw_q1,
      raw_q2 = raw_q2,
      local_varying = local_varying,
      total_variation = total_variation
    )
  )
}


.trsc_sim_weighted_static_simplex <- function(X, y, observation_weights,
                                              solver_control = list()) {
  if (!is.matrix(X) || !is.numeric(X) || !is.numeric(y) ||
      nrow(X) != length(y)) {
    stop("X and y must be aligned numeric fitting data.", call. = FALSE)
  }
  if (!is.numeric(observation_weights) ||
      length(observation_weights) != nrow(X) ||
      any(!is.finite(observation_weights)) ||
      any(observation_weights < 0) ||
      sum(observation_weights) <= 0) {
    stop("observation_weights must be finite nonnegative values with positive sum.",
         call. = FALSE)
  }
  J <- ncol(X)
  root_w <- sqrt(observation_weights)
  Xw <- X * root_w
  yw <- y * root_w

  P <- 2 * Matrix::Matrix(crossprod(Xw), sparse = TRUE)
  P <- Matrix::forceSymmetric(P, uplo = "U")
  qvec <- -2 * as.numeric(crossprod(Xw, yw))
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
    stop("Weighted simplex benchmark failed with OSQP status: ",
         result$info$status, call. = FALSE)
  }
  as.numeric(result$x)
}

trsc_sim_local_varying_simplex <- function(prepared, fitting_periods,
                                           bandwidth,
                                           solver_control = list()) {
  if (!inherits(prepared, "trsc_data")) {
    stop("prepared must be a trsc_data object.", call. = FALSE)
  }
  if (length(bandwidth) != 1L || !is.numeric(bandwidth) ||
      !is.finite(bandwidth) || bandwidth <= 0) {
    stop("bandwidth must be one finite strictly positive number.",
         call. = FALSE)
  }
  index <- match(fitting_periods, prepared$blocks$periods)
  if (anyNA(index)) stop("fitting_periods are not available.", call. = FALSE)
  X <- prepared$blocks$Y0[index, , drop = FALSE]
  y <- as.numeric(prepared$blocks$Y1[index])
  n <- nrow(X)
  J <- ncol(X)
  W <- matrix(NA_real_, nrow = J, ncol = n,
              dimnames = list(colnames(X), as.character(fitting_periods)))

  for (i in seq_len(n)) {
    distance <- as.numeric(fitting_periods - fitting_periods[i])
    kernel <- exp(-0.5 * (distance / bandwidth)^2)
    W[, i] <- .trsc_sim_weighted_static_simplex(
      X, y, kernel, solver_control = solver_control
    )
  }
  fitted <- rowSums(X * t(W))
  structure(
    list(
      weights = W,
      fitted = fitted,
      residuals = y - fitted,
      periods = fitting_periods,
      bandwidth = bandwidth,
      information = "same fitting outcomes as TRSC",
      benchmark = "local Gaussian-kernel varying-coefficient simplex smoother"
    ),
    class = "trsc_sim_local_varying_fit"
  )
}

trsc_sim_total_variation_fit <- function(prepared, fitting_periods, lambda,
                                         solver_control = list()) {
  if (!inherits(prepared, "trsc_data")) {
    stop("prepared must be a trsc_data object.", call. = FALSE)
  }
  if (length(lambda) != 1L || !is.numeric(lambda) ||
      !is.finite(lambda) || lambda < 0) {
    stop("lambda must be one finite nonnegative number.", call. = FALSE)
  }
  index <- match(fitting_periods, prepared$blocks$periods)
  if (anyNA(index)) stop("fitting_periods are not available.", call. = FALSE)
  X <- prepared$blocks$Y0[index, , drop = FALSE]
  y <- as.numeric(prepared$blocks$Y1[index])
  n <- nrow(X)
  J <- ncol(X)

  fit_blocks <- lapply(seq_len(n), function(i) {
    Matrix::Matrix(tcrossprod(X[i, ]), sparse = TRUE)
  })
  H_fit <- do.call(Matrix::bdiag, fit_blocks)
  m_w <- n * J
  m_d <- (n - 1L) * J
  P <- Matrix::bdiag(
    2 * H_fit,
    Matrix::Matrix(0, nrow = m_d, ncol = m_d, sparse = TRUE)
  )
  P <- Matrix::forceSymmetric(P, uplo = "U")
  qvec <- c(-2 * as.vector(t(X * y)), rep(lambda, m_d))

  D <- matrix(0, nrow = n - 1L, ncol = n)
  for (i in seq_len(n - 1L)) D[i, i:(i + 1L)] <- c(-1, 1)
  K <- kronecker(Matrix::Matrix(D, sparse = TRUE), Matrix::Diagonal(J))

  sum_row <- Matrix::Matrix(matrix(1, nrow = 1L, ncol = J), sparse = TRUE)
  simplex <- kronecker(Matrix::Diagonal(n), sum_row)
  zero_eq <- Matrix::Matrix(0, nrow = n, ncol = m_d, sparse = TRUE)
  A_eq <- cbind(simplex, zero_eq)

  Iw <- Matrix::Diagonal(m_w)
  Zu <- Matrix::Matrix(0, nrow = m_w, ncol = m_d, sparse = TRUE)
  A_w <- cbind(Iw, Zu)

  Zw <- Matrix::Matrix(0, nrow = m_d, ncol = m_w, sparse = TRUE)
  Iu <- Matrix::Diagonal(m_d)
  A_u <- cbind(Zw, Iu)

  A_pos <- cbind(-K, Iu)
  A_neg <- cbind(K, Iu)

  A <- rbind(A_eq, A_w, A_u, A_pos, A_neg)
  l <- c(
    rep(1, n),
    rep(0, m_w),
    rep(0, m_d),
    rep(0, m_d),
    rep(0, m_d)
  )
  u <- c(
    rep(1, n),
    rep(1, m_w),
    rep(Inf, m_d),
    rep(Inf, m_d),
    rep(Inf, m_d)
  )

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
    stop("Total-variation benchmark failed with OSQP status: ",
         result$info$status, call. = FALSE)
  }

  W <- matrix(result$x[seq_len(m_w)], nrow = J, ncol = n,
              dimnames = list(colnames(X), as.character(fitting_periods)))
  fitted <- rowSums(X * t(W))
  structure(
    list(
      weights = W,
      fitted = fitted,
      residuals = y - fitted,
      periods = fitting_periods,
      lambda = lambda,
      solver_status = result$info$status,
      information = "same fitting outcomes as TRSC",
      benchmark = "coordinate L1 total-variation simplex path"
    ),
    class = "trsc_sim_tv_fit"
  )
}

trsc_sim_historical_select <- function(prepared, candidates,
                                       validation_periods) {
  if (!inherits(prepared, "trsc_data")) {
    stop("prepared must be a trsc_data object.", call. = FALSE)
  }
  if (!is.list(candidates) || !length(candidates) ||
      !all(vapply(candidates, inherits, logical(1), what = "trsc_readout"))) {
    stop("candidates must be a nonempty frozen readout library.",
         call. = FALSE)
  }
  index <- match(validation_periods, prepared$blocks$periods)
  if (anyNA(index)) stop("validation_periods are not available.", call. = FALSE)
  donor_names <- colnames(prepared$blocks$Y0)
  Y0 <- prepared$blocks$Y0[index, , drop = FALSE]
  Y1 <- as.numeric(prepared$blocks$Y1[index])
  scores <- vapply(candidates, function(candidate) {
    w <- candidate$representative[donor_names]
    mean((Y1 - drop(Y0 %*% w))^2)
  }, numeric(1))
  selected_index <- which(scores == min(scores))[1L]
  structure(
    list(
      scores = scores,
      selected_index = selected_index,
      selected_name = names(candidates)[selected_index],
      selected = candidates[[selected_index]],
      target = "historical average validation risk",
      cutoff_targeted = FALSE,
      refit_performed = FALSE
    ),
    class = "trsc_sim_historical_selection"
  )
}

trsc_sim_design_diagnostics <- function(prepared, scenario, plan) {
  fit_index <- match(plan$fitting$periods, prepared$blocks$periods)
  X <- prepared$blocks$Y0[fit_index, , drop = FALSE]
  J <- ncol(X)
  gram <- crossprod(X) / nrow(X)
  gram_eig <- eigen((gram + t(gram)) / 2, symmetric = TRUE,
                    only.values = TRUE)$values

  contrast <- if (J > 1L) {
    sweep(X[, seq_len(J - 1L), drop = FALSE], 1L, X[, J], FUN = "-")
  } else {
    matrix(numeric(0), nrow = nrow(X), ncol = 0L)
  }
  contrast_eig <- if (ncol(contrast)) {
    eigen(crossprod(contrast) / nrow(contrast), symmetric = TRUE,
          only.values = TRUE)$values
  } else {
    numeric(0)
  }
  calibration_geometry <- trsc_calibration_geometry(prepared, plan)

  data.frame(
    scenario = scenario$scenario,
    fitting_dates = nrow(X),
    donors = J,
    fitting_gram_rank = qr(gram)$rank,
    fitting_gram_min_eigenvalue = min(gram_eig),
    contrast_rank = if (length(contrast_eig)) qr(contrast)$rank else 0L,
    contrast_min_eigenvalue = if (length(contrast_eig)) {
      min(contrast_eig)
    } else {
      NA_real_
    },
    calibration_rank = calibration_geometry$diagnostics$empirical_rank,
    reference_rank = qr(scenario$Gamma_reference)$rank,
    reference_min_eigenvalue = min(eigen(
      (scenario$Gamma_reference + t(scenario$Gamma_reference)) / 2,
      symmetric = TRUE, only.values = TRUE
    )$values),
    stringsAsFactors = FALSE
  )
}

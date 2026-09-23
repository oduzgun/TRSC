# Research-only simulation generators for TRSC.
#
# These functions are intentionally not part of the package namespace. They
# implement stress designs proposed by the methodological paper and return
# known data-generating objects for reproducible Monte Carlo work.

.trsc_sim_draw_gaussian <- function(n, Sigma) {
  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) {
    stop("Sigma must be a square matrix.", call. = FALSE)
  }
  eig <- eigen((Sigma + t(Sigma)) / 2, symmetric = TRUE)
  if (min(eig$values) < -1e-10) {
    stop("Sigma must be positive semidefinite.", call. = FALSE)
  }
  values <- pmax(eig$values, 0)
  root <- eig$vectors %*% diag(sqrt(values), nrow = length(values))
  matrix(stats::rnorm(n * nrow(Sigma)), nrow = n) %*% t(root)
}

.trsc_sim_path <- function(T, type) {
  u <- if (T == 1L) 0 else (seq_len(T) - 1) / (T - 1)
  switch(
    type,
    constant = rep(0.35, T),
    affine = 0.15 + 0.70 * u,
    curved = 0.50 + 0.30 * sin(2 * pi * u),
    abrupt = ifelse(u < 0.55, 0.20, 0.80),
    stop("Unknown path type.", call. = FALSE)
  )
}

.trsc_sim_default_plan_indices <- function(T) {
  if (T < 60L) stop("T must be at least 60 for the default simulation plan.", call. = FALSE)
  calibration_end <- floor(0.25 * T)
  fitting_start <- calibration_end + 3L
  fitting_end <- floor(0.70 * T)
  validation_start <- fitting_end + 3L
  validation_end <- T - 1L
  list(
    calibration = seq_len(calibration_end),
    fitting = fitting_start:fitting_end,
    validation = validation_start:validation_end,
    target = T
  )
}

#' Generate one research-only TRSC simulation scenario.
#'
#' @param scenario One of the paper-motivated design labels.
#' @param T Number of untreated dates.
#' @param J Number of donors before any exact duplicate is added.
#' @param seed Integer random seed.
#' @param noise_sd Standard deviation of target innovation noise.
#' @param reference_size Number of independent reference states used only for
#'   Monte Carlo evaluation of prediction/profile loss.
#'
#' @return A list containing a prepared-panel-ready data frame, true weight
#'   path, reference evaluation states and scenario metadata.
trsc_sim_scenario <- function(
    scenario = c(
      "constant", "affine", "curved", "abrupt",
      "exact_duplication", "near_redundancy",
      "weak_excitation", "reference_mismatch"
    ),
    T = 90L,
    J = 6L,
    seed = 1L,
    noise_sd = 0.20,
    reference_size = 4000L) {

  scenario <- match.arg(scenario)
  T <- as.integer(T)
  J <- as.integer(J)
  reference_size <- as.integer(reference_size)
  if (T < 60L || J < 3L || reference_size < 100L) {
    stop("Use T >= 60, J >= 3 and reference_size >= 100.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be one finite number.", call. = FALSE)
  }
  if (!is.numeric(noise_sd) || length(noise_sd) != 1L ||
      !is.finite(noise_sd) || noise_sd < 0) {
    stop("noise_sd must be one finite nonnegative number.", call. = FALSE)
  }
  set.seed(as.integer(seed))

  path_type <- if (scenario %in% c(
    "exact_duplication", "near_redundancy",
    "weak_excitation", "reference_mismatch"
  )) "curved" else scenario

  Sigma_reference <- diag(J)
  donor_names <- paste0("D", seq_len(J))

  reference_states <- .trsc_sim_draw_gaussian(reference_size, Sigma_reference)
  colnames(reference_states) <- donor_names

  donor_states <- .trsc_sim_draw_gaussian(T, Sigma_reference)
  colnames(donor_states) <- donor_names

  if (scenario == "weak_excitation") {
    plan_index <- .trsc_sim_default_plan_indices(T)
    common <- stats::rnorm(length(plan_index$fitting))
    weak <- outer(common, rep(1, J)) +
      matrix(stats::rnorm(length(plan_index$fitting) * J, sd = 0.03),
             ncol = J)
    donor_states[plan_index$fitting, ] <- weak
  }

  if (scenario == "reference_mismatch") {
    plan_index <- .trsc_sim_default_plan_indices(T)
    # Calibration/reference segment follows P. Fitting and validation donor
    # states are deliberately shifted/scaled to diagnose support mismatch.
    shifted <- c(plan_index$fitting, plan_index$validation, plan_index$target)
    donor_states[shifted, ] <- 1.25 * donor_states[shifted, , drop = FALSE]
    donor_states[shifted, 1] <- donor_states[shifted, 1] + 1.5
  }

  alpha <- .trsc_sim_path(T, path_type)
  W <- matrix(0, nrow = J, ncol = T,
              dimnames = list(donor_names, as.character(seq_len(T))))
  W[1, ] <- alpha
  W[2, ] <- 1 - alpha

  representation_identified <- TRUE

  if (scenario == "exact_duplication") {
    donor_states <- cbind(donor_states, D_duplicate = donor_states[, 1])
    reference_states <- cbind(reference_states,
                              D_duplicate = reference_states[, 1])
    W <- rbind(W, D_duplicate = rep(0, T))
    donor_names <- rownames(W)
    representation_identified <- FALSE
  }

  if (scenario == "near_redundancy") {
    donor_states[, J] <- donor_states[, 1] +
      stats::rnorm(T, sd = 0.02)
    reference_states[, J] <- reference_states[, 1] +
      stats::rnorm(reference_size, sd = 0.02)
  }

  target_mean <- rowSums(donor_states * t(W))
  noise <- stats::rnorm(T, sd = noise_sd)
  target <- target_mean + noise

  units <- c("Target", donor_names)
  panel <- expand.grid(
    unit = units,
    period = seq_len(T),
    stringsAsFactors = FALSE
  )
  lookup <- cbind(Target = target, donor_states)
  panel$outcome <- unlist(lapply(seq_len(T), function(t) {
    as.numeric(lookup[t, units])
  }), use.names = FALSE)

  plan_index <- .trsc_sim_default_plan_indices(T)
  Gamma_reference <- crossprod(reference_states) / nrow(reference_states)

  structure(
    list(
      scenario = scenario,
      seed = as.integer(seed),
      T = T,
      J = length(donor_names),
      donor_names = donor_names,
      panel = panel,
      donor_states = donor_states,
      target_mean = target_mean,
      target_outcome = target,
      true_weights = W,
      reference_states = reference_states,
      Gamma_reference = Gamma_reference,
      plan_index = plan_index,
      representation_identified = representation_identified,
      path_type = path_type,
      noise_sd = noise_sd,
      notes = switch(
        scenario,
        abrupt = "Abrupt path is intentionally outside the squared-difference smoothness class used for positive theory.",
        exact_duplication = "Raw donor representation is intentionally nonunique.",
        near_redundancy = "Last donor is a noisy near-copy of D1.",
        weak_excitation = "Fitting donor states are nearly collinear while the reference law remains full-dimensional.",
        reference_mismatch = "Calibration reference law differs from later fitting/validation donor-state support.",
        "Baseline smooth design."
      )
    ),
    class = "trsc_sim_scenario"
  )
}

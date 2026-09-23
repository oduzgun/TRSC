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



.trsc_sim_signs <- function(n, d) {
  matrix(sample(c(-1, 1), n * d, replace = TRUE), nrow = n, ncol = d)
}

.trsc_sim_face_vertices <- function(r, theta) {
  if (r < 3L) stop("r must be at least 3.", call. = FALSE)
  vertices <- matrix(0, nrow = 2L + 2L * (r - 2L), ncol = r)
  vertices[1L, 1L] <- -1
  vertices[2L, 1L] <- 1
  row <- 2L
  for (j in seq_len(r - 2L)) {
    row <- row + 1L
    vertices[row, 2L] <- 1
    vertices[row, j + 2L] <- theta
    row <- row + 1L
    vertices[row, 2L] <- 1
    vertices[row, j + 2L] <- -theta
  }
  rownames(vertices) <- c(
    "D_minus_u", "D_plus_u",
    unlist(lapply(seq_len(r - 2L), function(j) {
      c(paste0("D_v_plus_h", j), paste0("D_v_minus_h", j))
    }))
  )
  colnames(vertices) <- c("u", "v", paste0("h", seq_len(r - 2L)))
  vertices
}

.trsc_sim_dependent_intrinsic_states <- function(T, r, refresh_prob) {
  U <- sample(c(-1, 1), T, replace = TRUE)
  H <- matrix(0, nrow = T, ncol = r - 2L)
  H[1L, ] <- sample(c(-1, 1), r - 2L, replace = TRUE)
  if (T > 1L) {
    for (t in 2:T) {
      if (stats::runif(1) < refresh_prob) {
        H[t, ] <- sample(c(-1, 1), r - 2L, replace = TRUE)
      } else {
        H[t, ] <- H[t - 1L, ]
      }
    }
  }
  cbind(u = U, v = 1, H)
}

.trsc_sim_face_scenario <- function(scenario, T, J, seed, noise_sd,
                                    reference_size,
                                    theta = 1 / 16,
                                    face_margin = NULL,
                                    refresh_prob = 0.35) {
  set.seed(as.integer(seed))
  r <- max(3L, as.integer(ceiling((J + 2L) / 2L)))
  vertices <- .trsc_sim_face_vertices(r, theta)
  donor_names <- rownames(vertices)
  J_actual <- nrow(vertices)

  donor_z <- .trsc_sim_dependent_intrinsic_states(T, r, refresh_prob)
  donor_states <- donor_z %*% t(vertices)
  colnames(donor_states) <- donor_names

  reference_z <- cbind(
    u = sample(c(-1, 1), reference_size, replace = TRUE),
    v = 1,
    .trsc_sim_signs(reference_size, r - 2L)
  )
  heldout_z <- cbind(
    u = sample(c(-1, 1), reference_size, replace = TRUE),
    v = 1,
    .trsc_sim_signs(reference_size, r - 2L)
  )
  colnames(reference_z) <- colnames(vertices)
  colnames(heldout_z) <- colnames(vertices)
  reference_states <- reference_z %*% t(vertices)
  heldout_states <- heldout_z %*% t(vertices)
  colnames(reference_states) <- donor_names
  colnames(heldout_states) <- donor_names

  b <- 0.55 * sin(2 * pi * (seq_len(T) - 1) / max(1, T - 1))
  margin <- if (scenario == "vanishing_margin") {
    if (is.null(face_margin)) 0.02 else face_margin
  } else {
    if (is.null(face_margin)) 0.40 else face_margin
  }
  if (!is.finite(margin) || margin <= 0) {
    stop("face_margin must be strictly positive.", call. = FALSE)
  }

  eps_bound <- min(margin / 4, sqrt(3) * noise_sd)
  epsilon <- if (eps_bound > 0) {
    stats::runif(T, min = -eps_bound, max = eps_bound)
  } else {
    rep(0, T)
  }

  W <- matrix(0, nrow = J_actual, ncol = T,
              dimnames = list(donor_names, as.character(seq_len(T))))
  W[1L, ] <- (1 - b) / 2
  W[2L, ] <- (1 + b) / 2

  target_mean <- b * donor_z[, "u"] - margin
  target <- target_mean + epsilon
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
  Gamma_heldout <- crossprod(heldout_states) / nrow(heldout_states)

  structure(
    list(
      scenario = scenario,
      seed = as.integer(seed),
      T = T,
      J = J_actual,
      donor_names = donor_names,
      panel = panel,
      donor_states = donor_states,
      target_mean = target_mean,
      target_outcome = target,
      true_weights = W,
      reference_states = reference_states,
      heldout_states = heldout_states,
      Gamma_reference = Gamma_reference,
      Gamma_heldout = Gamma_heldout,
      plan_index = plan_index,
      representation_identified = TRUE,
      path_type = "smooth_separated_edge",
      noise_sd = stats::sd(epsilon),
      structural = list(
        intrinsic_dimension = r,
        theta = theta,
        face_margin = margin,
        refresh_prob = refresh_prob,
        target_face = "edge{-u,+u}",
        positive_normal_margin = scenario == "separated_face",
        theorem_class = scenario == "separated_face"
      ),
      notes = if (scenario == "separated_face") {
        "Paper-matched stable exposed-edge stress design with positive outward normal margin."
      } else {
        "Same exposed-edge geometry with a deliberately small normal margin; used to stress the face-localization boundary."
      }
    ),
    class = "trsc_sim_scenario"
  )
}

.trsc_sim_face_change_scenario <- function(T, seed, noise_sd, reference_size) {
  set.seed(as.integer(seed))
  donor_names <- c("D_bottom_left", "D_bottom_right",
                   "D_top_left", "D_top_right")
  vertices <- rbind(
    D_bottom_left = c(-1, 0),
    D_bottom_right = c(1, 0),
    D_top_left = c(-1, 1),
    D_top_right = c(1, 1)
  )
  reference_z <- .trsc_sim_draw_gaussian(reference_size, diag(2))
  heldout_z <- .trsc_sim_draw_gaussian(reference_size, diag(2))
  donor_z <- .trsc_sim_draw_gaussian(T, diag(2))
  donor_states <- donor_z %*% t(vertices)
  reference_states <- reference_z %*% t(vertices)
  heldout_states <- heldout_z %*% t(vertices)
  colnames(donor_states) <- donor_names
  colnames(reference_states) <- donor_names
  colnames(heldout_states) <- donor_names

  W <- matrix(0, nrow = 4L, ncol = T,
              dimnames = list(donor_names, as.character(seq_len(T))))
  switch_t <- floor(T * 0.58)
  b <- 0.45 * sin(2 * pi * (seq_len(T) - 1) / max(1, T - 1))
  bottom_left <- (1 - b) / 2
  bottom_right <- (1 + b) / 2
  W[1L, seq_len(switch_t)] <- bottom_left[seq_len(switch_t)]
  W[2L, seq_len(switch_t)] <- bottom_right[seq_len(switch_t)]
  after <- (switch_t + 1L):T
  top_left <- (1 - b[after]) / 2
  top_right <- (1 + b[after]) / 2
  W[3L, after] <- top_left
  W[4L, after] <- top_right

  target_mean <- rowSums(donor_states * t(W))
  target <- target_mean + stats::rnorm(T, sd = noise_sd)
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

  structure(
    list(
      scenario = "face_change",
      seed = as.integer(seed),
      T = T,
      J = 4L,
      donor_names = donor_names,
      panel = panel,
      donor_states = donor_states,
      target_mean = target_mean,
      target_outcome = target,
      true_weights = W,
      reference_states = reference_states,
      heldout_states = heldout_states,
      Gamma_reference = crossprod(reference_states) / nrow(reference_states),
      Gamma_heldout = crossprod(heldout_states) / nrow(heldout_states),
      plan_index = .trsc_sim_default_plan_indices(T),
      representation_identified = TRUE,
      path_type = "changing_faces",
      noise_sd = noise_sd,
      structural = list(
        switch_period = switch_t,
        face_before = "bottom edge",
        face_after = "top edge",
        theorem_class = FALSE
      ),
      notes = "Exact-representation stress design that switches between two exposed edges; no generic face-change theorem is claimed."
    ),
    class = "trsc_sim_scenario"
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
      "weak_excitation", "reference_mismatch",
      "interior", "persistent_errors",
      "separated_face", "vanishing_margin", "face_change"
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
  if (scenario %in% c("separated_face", "vanishing_margin")) {
    return(.trsc_sim_face_scenario(
      scenario = scenario, T = T, J = J, seed = seed,
      noise_sd = noise_sd, reference_size = reference_size
    ))
  }
  if (scenario == "face_change") {
    return(.trsc_sim_face_change_scenario(
      T = T, seed = seed, noise_sd = noise_sd,
      reference_size = reference_size
    ))
  }
  set.seed(as.integer(seed))

  path_type <- if (scenario %in% c(
    "exact_duplication", "near_redundancy",
    "weak_excitation", "reference_mismatch",
    "interior", "persistent_errors"
  )) "curved" else scenario

  Sigma_reference <- diag(J)
  donor_names <- paste0("D", seq_len(J))

  reference_states <- .trsc_sim_draw_gaussian(reference_size, Sigma_reference)
  colnames(reference_states) <- donor_names

  heldout_states <- .trsc_sim_draw_gaussian(reference_size, Sigma_reference)
  colnames(heldout_states) <- donor_names

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
    # states are deliberately shifted/scaled. Independent held-out states
    # follow the same shifted law Q so evaluation can distinguish reference-P
    # profile loss from deployment/support-shift prediction loss.
    shifted <- c(plan_index$fitting, plan_index$validation, plan_index$target)
    donor_states[shifted, ] <- 1.25 * donor_states[shifted, , drop = FALSE]
    donor_states[shifted, 1] <- donor_states[shifted, 1] + 1.5
    heldout_states <- 1.25 * heldout_states
    heldout_states[, 1] <- heldout_states[, 1] + 1.5
  }

  alpha <- .trsc_sim_path(T, path_type)
  W <- matrix(0, nrow = J, ncol = T,
              dimnames = list(donor_names, as.character(seq_len(T))))
  if (scenario == "interior") {
    W[,] <- 1 / J
    interior_shift <- (0.30 / J) *
      sin(2 * pi * (seq_len(T) - 1) / max(1, T - 1))
    W[1, ] <- W[1, ] + interior_shift
    W[2, ] <- W[2, ] - interior_shift
  } else {
    W[1, ] <- alpha
    W[2, ] <- 1 - alpha
  }

  representation_identified <- TRUE

  if (scenario == "exact_duplication") {
    donor_states <- cbind(donor_states, D_duplicate = donor_states[, 1])
    reference_states <- cbind(reference_states,
                              D_duplicate = reference_states[, 1])
    heldout_states <- cbind(heldout_states,
                            D_duplicate = heldout_states[, 1])
    W <- rbind(W, D_duplicate = rep(0, T))
    donor_names <- rownames(W)
    representation_identified <- FALSE
  }

  if (scenario == "near_redundancy") {
    donor_states[, J] <- donor_states[, 1] +
      stats::rnorm(T, sd = 0.02)
    reference_states[, J] <- reference_states[, 1] +
      stats::rnorm(reference_size, sd = 0.02)
    heldout_states[, J] <- heldout_states[, 1] +
      stats::rnorm(reference_size, sd = 0.02)
  }

  target_mean <- rowSums(donor_states * t(W))
  if (scenario == "persistent_errors") {
    innovations <- stats::rnorm(T, sd = noise_sd)
    noise <- numeric(T)
    rho <- 0.85
    noise[1] <- innovations[1] / sqrt(max(1e-8, 1 - rho^2))
    if (T > 1L) {
      for (t in 2:T) noise[t] <- rho * noise[t - 1L] + innovations[t]
    }
  } else {
    noise <- stats::rnorm(T, sd = noise_sd)
  }
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
  Gamma_heldout <- crossprod(heldout_states) / nrow(heldout_states)

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
      heldout_states = heldout_states,
      Gamma_reference = Gamma_reference,
      Gamma_heldout = Gamma_heldout,
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
        reference_mismatch = "Calibration reference law P differs from later fitting/validation and held-out deployment law Q.",
        interior = "Target path remains strictly inside the donor simplex; this stresses the full-dimensional tangent regime.",
        persistent_errors = "Target innovation follows a persistent AR(1) process, inflating temporal covariance outside the independent-noise benchmark.",
        "Baseline smooth design."
      )
    ),
    class = "trsc_sim_scenario"
  )
}

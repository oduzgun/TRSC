#' Compute the dependent-panel coupling bound for a TRSC sample plan
#'
#' @param plan A `trsc_sample_plan` object.
#' @param beta A function returning the assumed absolute-regularity coefficient
#'   `beta(k)` for a positive integer separation `k`.
#'
#' @return A `trsc_mixing_bound` object containing the four contributions and
#'   their sum.
#'
#' @details For retained calibration count `N_C`, retained validation count
#'   `v`, thinning spacings `b_C` and `b_V`, and effective gaps `g_0`
#'   and `g_1`, the reported quantity is
#'
#'   `(N_C - 1) beta(b_C) + beta(g_0) + beta(g_1) + (v - 1) beta(b_V)`.
#'
#'   This is a calculation under a user-supplied mixing-coefficient function.
#'   The package does not estimate or certify that function, stationarity, or
#'   the other fitting and validation assumptions.
#'
#' @export
trsc_mixing_bound <- function(plan, beta) {
  if (!inherits(plan, "trsc_sample_plan") || !is.list(plan)) {
    stop("plan must be a trsc_sample_plan object.", call. = FALSE)
  }
  if (!is.function(beta)) {
    stop("beta must be a function beta(k).", call. = FALSE)
  }

  bC <- plan$specification$calibration_thin
  bV <- plan$specification$validation_thin
  g0 <- plan$gaps$g0_steps
  g1 <- plan$gaps$g1_steps
  NC <- plan$calibration$retained_count
  v <- plan$validation$retained_count

  evaluate_beta <- function(k, label) {
    value <- beta(k)
    if (length(value) != 1L || !is.numeric(value) || is.complex(value) ||
        is.na(value) || !is.finite(value) || value < 0 || value > 1) {
      stop("beta(", label, ") must return one finite value in [0, 1].",
           call. = FALSE)
    }
    as.numeric(value)
  }

  beta_bC <- evaluate_beta(bC, "b_C")
  beta_bV <- evaluate_beta(bV, "b_V")
  beta_g0 <- evaluate_beta(g0, "g_0")
  beta_g1 <- evaluate_beta(g1, "g_1")

  contributions <- c(
    calibration_internal = (NC - 1L) * beta_bC,
    calibration_to_fitting = beta_g0,
    fitting_to_validation = beta_g1,
    validation_internal = (v - 1L) * beta_bV
  )

  structure(
    list(
      schema_version = 1L,
      retained_calibration = NC,
      retained_validation = v,
      calibration_thin = bC,
      validation_thin = bV,
      gap0 = g0,
      gap1 = g1,
      beta_values = c(
        beta_calibration_thin = beta_bC,
        beta_gap0 = beta_g0,
        beta_gap1 = beta_g1,
        beta_validation_thin = beta_bV
      ),
      contributions = contributions,
      delta_mix = sum(contributions),
      certification = list(
        beta_function_verified = FALSE,
        stationarity_verified = FALSE,
        absolute_regularity_verified = FALSE,
        fitting_assumptions_verified = FALSE,
        validation_assumptions_verified = FALSE
      )
    ),
    class = "trsc_mixing_bound"
  )
}

#' Construct the paper's explicit exponential-mixing spacing design
#'
#' @param calibration_raw_length Positive integer raw calibration-segment
#'   length.
#' @param validation_raw_length Positive integer raw validation-segment length.
#' @param C_beta Positive finite constant in the assumed envelope
#'   `beta(k) <= C_beta * exp(-k / tau_beta)`.
#' @param tau_beta Positive finite decay scale.
#' @param delta_mix Desired finite coupling-error budget in `(0, 1)`.
#'
#' @return A `trsc_geometric_mixing_design` containing the recommended
#'   calibration and validation thinning spacings, minimum gaps, implied
#'   retained counts, and the algebraic envelope bound.
#'
#' @details This helper implements the explicit conservative choices in the
#'   pre-treatment dependence-transfer proposition. It does not estimate
#'   `C_beta` or `tau_beta` and does not certify geometric beta-mixing.
#'   If the implied validation retention contains fewer than two dates, the
#'   requested raw validation segment is insufficient for cutoff validation.
#'
#' @export
trsc_geometric_mixing_design <- function(calibration_raw_length,
                                         validation_raw_length,
                                         C_beta,
                                         tau_beta,
                                         delta_mix) {
  calibration_raw_length <- .trsc_positive_integer(
    calibration_raw_length, "calibration_raw_length"
  )
  validation_raw_length <- .trsc_positive_integer(
    validation_raw_length, "validation_raw_length"
  )
  if (length(C_beta) != 1L || !is.numeric(C_beta) || is.complex(C_beta) ||
      is.na(C_beta) || !is.finite(C_beta) || C_beta <= 0) {
    stop("C_beta must be one finite strictly positive number.", call. = FALSE)
  }
  if (length(tau_beta) != 1L || !is.numeric(tau_beta) || is.complex(tau_beta) ||
      is.na(tau_beta) || !is.finite(tau_beta) || tau_beta <= 0) {
    stop("tau_beta must be one finite strictly positive number.", call. = FALSE)
  }
  if (length(delta_mix) != 1L || !is.numeric(delta_mix) ||
      is.complex(delta_mix) || is.na(delta_mix) ||
      !is.finite(delta_mix) || delta_mix <= 0 || delta_mix >= 1) {
    stop("delta_mix must lie strictly between 0 and 1.", call. = FALSE)
  }

  bC <- max(
    1L,
    as.integer(ceiling(
      tau_beta * log(4 * C_beta * calibration_raw_length / delta_mix)
    ))
  )
  bV <- max(
    1L,
    as.integer(ceiling(
      tau_beta * log(4 * C_beta * validation_raw_length / delta_mix)
    ))
  )
  gap <- max(
    1L,
    as.integer(ceiling(
      tau_beta * log(4 * C_beta / delta_mix)
    ))
  )

  NC <- 1L + floor((calibration_raw_length - 1L) / bC)
  v <- 1L + floor((validation_raw_length - 1L) / bV)
  if (v < 2L) {
    stop(
      "The requested validation raw length retains fewer than two dates under ",
      "the explicit geometric-mixing spacing. Increase validation_raw_length ",
      "or relax the supplied mixing/error parameters.",
      call. = FALSE
    )
  }

  envelope <- function(k) C_beta * exp(-k / tau_beta)
  contributions <- c(
    calibration_internal = (NC - 1L) * envelope(bC),
    calibration_to_fitting = envelope(gap),
    fitting_to_validation = envelope(gap),
    validation_internal = (v - 1L) * envelope(bV)
  )
  upper_bound <- sum(contributions)

  structure(
    list(
      schema_version = 1L,
      assumption = list(
        form = "beta(k) <= C_beta * exp(-k / tau_beta)",
        C_beta = C_beta,
        tau_beta = tau_beta,
        certified = FALSE
      ),
      delta_mix_target = delta_mix,
      calibration_raw_length = calibration_raw_length,
      validation_raw_length = validation_raw_length,
      calibration_thin = bC,
      validation_thin = bV,
      min_gap0 = gap,
      min_gap1 = gap,
      retained_calibration = as.integer(NC),
      retained_validation = as.integer(v),
      envelope_contributions = contributions,
      envelope_upper_bound = upper_bound,
      target_respected_numerically =
        upper_bound <= delta_mix * (1 + 100 * .Machine$double.eps)
    ),
    class = "trsc_geometric_mixing_design"
  )
}

#' @export
print.trsc_mixing_bound <- function(x, ...) {
  cat("TRSC dependent-panel coupling diagnostic\n")
  cat("  Retained calibration: ", x$retained_calibration, "\n", sep = "")
  cat("  Retained validation: ", x$retained_validation, "\n", sep = "")
  cat("  Delta_mix: ", format(x$delta_mix), "\n", sep = "")
  cat("  Mixing law certified by software: no\n")
  invisible(x)
}

#' @export
print.trsc_geometric_mixing_design <- function(x, ...) {
  cat("TRSC explicit geometric-mixing spacing design\n")
  cat("  Calibration thinning: ", x$calibration_thin, "\n", sep = "")
  cat("  Validation thinning: ", x$validation_thin, "\n", sep = "")
  cat("  Minimum gaps g0 / g1: ", x$min_gap0, " / ", x$min_gap1,
      "\n", sep = "")
  cat("  Retained calibration / validation: ", x$retained_calibration,
      " / ", x$retained_validation, "\n", sep = "")
  cat("  Envelope upper bound: ", format(x$envelope_upper_bound), "\n", sep = "")
  cat("  Mixing assumption certified by software: no\n")
  invisible(x)
}

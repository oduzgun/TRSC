#' Construct a frozen terminal readout from a TRSC fit
#'
#' @param fit A `trsc_fit` object whose fitted periods end at the fitting
#'   cutoff.
#' @param h Positive integer terminal bandwidth, at most the number of fitted
#'   dates.
#'
#' @return A `trsc_readout` object. The readout is a convex taper of the last
#'   `h` fitted relationships and is frozen once constructed.
#'
#' @details For lags `j = 0, ..., h - 1`, the taper is
#'   `2 * (h - j) / (h * (h + 1))`. Thus `h = 1` is the endpoint readout.
#'   The same convex average is applied to one numerical simplex lift and to
#'   supplied-geometry metric coordinates. The numerical lift remains a
#'   computational representative rather than an identified structural share.
#'
#' @export
trsc_readout <- function(fit, h = 1L) {
  if (!inherits(fit, "trsc_fit") || !is.list(fit)) {
    stop("fit must be a trsc_fit object.", call. = FALSE)
  }
  if (length(h) != 1L || !is.numeric(h) || is.complex(h) ||
      is.na(h) || !is.finite(h) || h != as.integer(h) || h < 1L) {
    stop("h must be one positive integer-valued number.", call. = FALSE)
  }
  h <- as.integer(h)
  n <- length(fit$periods)
  if (h > n) {
    stop("h cannot exceed the number of fitted dates.", call. = FALSE)
  }

  lags <- 0:(h - 1L)
  taper <- 2 * (h - lags) / (h * (h + 1))
  terminal_index <- n - lags

  W <- fit$weights[, terminal_index, drop = FALSE]
  representative <- drop(W %*% taper)
  names(representative) <- fit$donor_names

  metric <- fit$metric_coordinates[, terminal_index, drop = FALSE]
  metric_coordinates <- if (nrow(metric)) {
    drop(metric %*% taper)
  } else {
    numeric(0)
  }
  if (length(metric_coordinates)) {
    names(metric_coordinates) <- rownames(fit$metric_coordinates)
  }

  structure(
    list(
      schema_version = 1L,
      bandwidth = h,
      cutoff = fit$periods[n],
      source_periods = fit$periods,
      source_specification = fit$specification,
      donor_names = fit$donor_names,
      taper_weights = stats::setNames(taper, as.character(fit$periods[terminal_index])),
      representative = representative,
      metric_coordinates = metric_coordinates,
      geometry = fit$geometry,
      frozen = TRUE
    ),
    class = "trsc_readout"
  )
}

.trsc_default_bandwidths <- function(n) {
  if (length(n) != 1L || !is.numeric(n) || is.complex(n) ||
      is.na(n) || !is.finite(n) || n != as.integer(n) || n < 1L) {
    stop("n must be one positive integer-valued number.", call. = FALSE)
  }
  n <- as.integer(n)
  powers <- 2^(0:floor(log(n, base = 2)))
  as.integer(sort(unique(c(powers[powers <= n], n))))
}

#' Build the default finite terminal-readout library
#'
#' @param fit A `trsc_fit` object.
#' @param bandwidths Optional vector of distinct positive integer bandwidths.
#'   If NULL, use `1, 2, 4, ...` together with the full fitted length.
#'
#' @return A named list of frozen `trsc_readout` objects with class
#'   `trsc_readout_library`.
#' @export
trsc_readout_library <- function(fit, bandwidths = NULL) {
  if (!inherits(fit, "trsc_fit")) {
    stop("fit must be a trsc_fit object.", call. = FALSE)
  }
  n <- length(fit$periods)
  if (is.null(bandwidths)) {
    bandwidths <- .trsc_default_bandwidths(n)
  } else {
    if (!is.numeric(bandwidths) || is.complex(bandwidths) ||
        !length(bandwidths) || anyNA(bandwidths) ||
        any(!is.finite(bandwidths)) ||
        any(bandwidths != as.integer(bandwidths)) ||
        any(bandwidths < 1L) || any(bandwidths > n) ||
        anyDuplicated(bandwidths)) {
      stop("bandwidths must be distinct positive integers not exceeding the fitted length.",
           call. = FALSE)
    }
    bandwidths <- as.integer(bandwidths)
  }
  out <- lapply(bandwidths, function(h) trsc_readout(fit, h))
  names(out) <- paste0("h", bandwidths)
  class(out) <- c("trsc_readout_library", "list")
  out
}

.trsc_local_linear_validation_weights <- function(periods, target_period,
                                                   kernel_weights = NULL,
                                                   tolerance = 1e-12) {
  if (!is.numeric(periods) || is.complex(periods) || length(periods) < 2L ||
      anyNA(periods) || any(!is.finite(periods)) || anyDuplicated(periods) ||
      any(diff(periods) <= 0)) {
    stop("validation_periods must contain at least two strictly increasing finite dates.",
         call. = FALSE)
  }
  if (length(target_period) != 1L || !is.numeric(target_period) ||
      is.complex(target_period) || is.na(target_period) ||
      !is.finite(target_period)) {
    stop("target_period must be one finite numeric date.", call. = FALSE)
  }
  if (is.null(kernel_weights)) {
    kernel_weights <- rep(1, length(periods))
  }
  if (!is.numeric(kernel_weights) || is.complex(kernel_weights) ||
      length(kernel_weights) != length(periods) || anyNA(kernel_weights) ||
      any(!is.finite(kernel_weights)) || any(kernel_weights < 0) ||
      sum(kernel_weights > 0) < 2L) {
    stop("kernel_weights must be finite nonnegative values with at least two positive entries.",
         call. = FALSE)
  }

  d <- periods - target_period
  S0 <- sum(kernel_weights)
  S1 <- sum(kernel_weights * d)
  S2 <- sum(kernel_weights * d^2)
  denominator <- S0 * S2 - S1^2
  scale <- max(1, abs(S0 * S2), S1^2)
  if (!is.finite(denominator) || denominator <= tolerance * scale) {
    stop("Validation dates and kernel weights do not identify a local-linear readout.",
         call. = FALSE)
  }

  omega <- kernel_weights * (S2 - S1 * d) / denominator
  if (abs(sum(omega) - 1) > 100 * tolerance ||
      abs(sum(omega * d)) > 100 * tolerance * max(1, max(abs(d)))) {
    stop("Numerical failure while constructing cutoff-targeted validation weights.",
         call. = FALSE)
  }
  stats::setNames(omega, as.character(periods))
}

#' Select among a frozen TRSC readout library at the untreated cutoff
#'
#' @param prepared A schema-version-2 `trsc_data` object containing the
#'   validation dates.
#' @param candidates A nonempty list of frozen `trsc_readout` objects.
#' @param validation_periods At least two strictly increasing pre-treatment
#'   dates occurring after every candidate's fitting cutoff.
#' @param target_period Untreated cutoff date. By default, and for the current
#'   theorem-aligned implementation, this must equal the last declared
#'   pre-treatment date.
#' @param kernel_weights Optional nonnegative deterministic kernel weights used
#'   to construct signed local-linear validation weights. The software cannot
#'   verify that user-supplied kernel choices were made without inspecting
#'   validation outcomes.
#'
#' @return A `trsc_selection` object containing the selected frozen readout,
#'   candidate scores, cutoff-targeted validation weights, and diagnostics.
#'
#' @details Candidates are scored by the signed weighted squared residual
#'   criterion from the frozen cutoff-selection procedure. This function does
#'   not refit candidates and does not perform a final full-sample refit.
#'   It checks chronological separation mechanically, but it does not certify
#'   the probabilistic independence, mixing, smoothness, or calibration
#'   assumptions used by the paper's oracle inequality.
#'
#' @export
trsc_select <- function(prepared, candidates, validation_periods,
                        target_period = NULL, kernel_weights = NULL) {
  if (!inherits(prepared, "trsc_data") || !is.list(prepared) ||
      !identical(prepared$schema_version, 2L) ||
      !inherits(prepared$blocks, "trsc_blocks")) {
    stop("prepared must be a schema-version-2 trsc_data object.", call. = FALSE)
  }
  if (!is.list(candidates) || !length(candidates) ||
      !all(vapply(candidates, inherits, logical(1), what = "trsc_readout"))) {
    stop("candidates must be a nonempty list of trsc_readout objects.",
         call. = FALSE)
  }

  candidate_names <- names(candidates)
  if (is.null(candidate_names)) {
    candidate_names <- sprintf("candidate_%02d", seq_along(candidates))
  }
  if (anyNA(candidate_names) || any(!nzchar(candidate_names)) ||
      anyDuplicated(candidate_names)) {
    stop("candidate names must be unique and nonempty.", call. = FALSE)
  }

  all_periods <- prepared$blocks$periods
  if (!is.numeric(validation_periods) || is.complex(validation_periods) ||
      length(validation_periods) < 2L || anyNA(validation_periods) ||
      any(!is.finite(validation_periods)) || anyDuplicated(validation_periods) ||
      any(diff(validation_periods) <= 0)) {
    stop("validation_periods must contain at least two strictly increasing finite dates.",
         call. = FALSE)
  }
  validation_index <- match(validation_periods, all_periods)
  if (anyNA(validation_index)) {
    stop("Every validation period must belong to the declared pre-treatment period.",
         call. = FALSE)
  }

  cutoff <- max(all_periods)
  if (is.null(target_period)) target_period <- cutoff
  if (length(target_period) != 1L || !is.numeric(target_period) ||
      is.complex(target_period) || is.na(target_period) ||
      !is.finite(target_period) || target_period != cutoff) {
    stop("target_period must equal the last declared pre-treatment date.",
         call. = FALSE)
  }
  if (max(validation_periods) > target_period) {
    stop("validation periods cannot occur after target_period.", call. = FALSE)
  }

  donor_names <- colnames(prepared$blocks$Y0)
  common_source_periods <- candidates[[1L]]$source_periods
  common_cutoff <- candidates[[1L]]$cutoff
  for (i in seq_along(candidates)) {
    candidate <- candidates[[i]]
    if (!identical(candidate$source_periods, common_source_periods) ||
        !identical(candidate$cutoff, common_cutoff)) {
      stop("All candidates must be frozen from the same fitting block and cutoff.",
           call. = FALSE)
    }
    if (!isTRUE(candidate$frozen)) {
      stop("Every candidate must be frozen before validation.", call. = FALSE)
    }
    if (candidate$cutoff >= min(validation_periods)) {
      stop("Every candidate fitting cutoff must precede all validation periods.",
           call. = FALSE)
    }
    if (!setequal(candidate$donor_names, donor_names)) {
      stop("Every candidate donor set must match prepared donors exactly.",
           call. = FALSE)
    }
    w <- candidate$representative
    if (!is.numeric(w) || any(!is.finite(w)) ||
        length(w) != length(donor_names) ||
        is.null(names(w)) || !setequal(names(w), donor_names)) {
      stop("Each candidate must contain a finite named donor representative.",
           call. = FALSE)
    }
  }

  omega <- .trsc_local_linear_validation_weights(
    validation_periods, target_period, kernel_weights
  )
  Y1 <- as.numeric(prepared$blocks$Y1[validation_index])
  Y0 <- prepared$blocks$Y0[validation_index, , drop = FALSE]

  scores <- vapply(seq_along(candidates), function(i) {
    w <- candidates[[i]]$representative[donor_names]
    prediction <- drop(Y0 %*% w)
    sum(omega * (Y1 - prediction)^2)
  }, numeric(1))
  names(scores) <- candidate_names

  minimum_score <- min(scores)
  tied_index <- which(scores == minimum_score)
  selected_index <- tied_index[1L]
  selected <- candidates[[selected_index]]
  names(candidates) <- candidate_names
  tie_rule <- if (inherits(candidates, "trsc_candidate_library") &&
                  !is.null(attr(candidates, "tie_rule"))) {
    attr(candidates, "tie_rule")
  } else {
    "first exact minimum in candidate order"
  }

  d <- validation_periods - target_period
  structure(
    list(
      schema_version = 1L,
      target_period = target_period,
      validation_periods = validation_periods,
      validation_weights = omega,
      kernel_weights = if (is.null(kernel_weights)) {
        stats::setNames(rep(1, length(validation_periods)),
                        as.character(validation_periods))
      } else {
        stats::setNames(as.numeric(kernel_weights),
                        as.character(validation_periods))
      },
      candidates = candidates,
      scores = scores,
      selected_index = selected_index,
      selected_name = candidate_names[selected_index],
      selected = selected,
      diagnostics = list(
        weight_sum = sum(omega),
        first_moment = sum(omega * d),
        signed_weights = any(omega < 0),
        chronological_separation_checked = TRUE,
        candidates_frozen = TRUE,
        exact_tie_count = length(tied_index),
        exact_tie_names = candidate_names[tied_index],
        tie_rule = tie_rule,
        refit_performed = FALSE,
        theorem_assumptions_certified = FALSE
      )
    ),
    class = "trsc_selection"
  )
}

#' @export
print.trsc_readout <- function(x, ...) {
  cat("Frozen TRSC terminal readout\n")
  cat("  Fitting cutoff: ", x$cutoff, "\n", sep = "")
  cat("  Bandwidth h: ", x$bandwidth, "\n", sep = "")
  cat("  Donors: ", length(x$donor_names), "\n", sep = "")
  cat("  Frozen before validation: ", if (isTRUE(x$frozen)) "yes" else "no",
      "\n", sep = "")
  invisible(x)
}

#' @export
print.trsc_selection <- function(x, ...) {
  cat("TRSC frozen cutoff selection\n")
  cat("  Target period: ", x$target_period, "\n", sep = "")
  cat("  Validation dates: ", length(x$validation_periods), "\n", sep = "")
  cat("  Candidates: ", length(x$candidates), "\n", sep = "")
  cat("  Selected: ", x$selected_name, "\n", sep = "")
  cat("  Refit performed: no\n")
  cat("  Theorem assumptions certified: no\n")
  invisible(x)
}

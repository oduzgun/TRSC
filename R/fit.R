#' Fit a pre-treatment Temporally Regularized Synthetic Control
#'
#' Estimate a path of simplex donor representations using squared fitting loss
#' plus first- or second-order temporal roughness measured by an explicitly
#' supplied predictive geometry.
#'
#' @param prepared A schema-version-2 `trsc_data` object.
#' @param geometry A `trsc_geometry` object constructed from a reference or
#'   calibration sample whose donor columns correspond exactly to the fitted
#'   donor set. The function does not assume that the fitting observations
#'   themselves define the reference law.
#' @param kappa One finite nonnegative regularization parameter.
#' @param q Temporal difference order, exactly 1 or 2.
#' @param periods Optional contiguous subset of declared pre-treatment periods.
#'   If NULL, all declared pre-treatment periods are used.
#' @param solver Currently only `"osqp"`.
#' @param solver_control Named list passed to `osqp::osqpSettings()`. Defaults
#'   use tight absolute and relative tolerances and polishing.
#' @param allow_inaccurate Logical. If FALSE, the OSQP status
#'   `"solved inaccurate"` is rejected.
#' @param feasibility_tolerance Finite nonnegative tolerance used only to
#'   verify simplex feasibility of the numerical solution. The returned
#'   weights are not clipped or renormalized.
#'
#' @return An object of class `trsc_fit`. Its primary fitted representation is
#'   the predictive path stored in empirical quotient coordinates. Raw donor
#'   weights are retained as one solver-selected simplex representation and
#'   need not be unique.
#'
#' @details The fitted objective is
#'   deqn{\sum_t (y_t-x_t^\top w_t)^2 +
#'   \kappa\sum_{t=q+1}^n (\Delta^q w_t)^\top\widehat\Gamma
#'   (\Delta^q w_t),}
#'   subject to eqn{w_t} lying in the donor simplex at every fitted date.
#'   The supplied geometry may be singular. No inverse of its second-moment
#'   matrix is formed and no coordinate ridge is added to force uniqueness.
#'
#'   Supplying a geometry estimated from the same fitting observations is
#'   computationally permitted, but this function does not certify that such a
#'   choice satisfies the sample-separation or calibration conditions of the
#'   statistical theory.
#'
#' @export
trsc_fit <- function(prepared, geometry, kappa, q = 1,
                     periods = NULL,
                     solver = c("osqp"),
                     solver_control = list(),
                     allow_inaccurate = FALSE,
                     feasibility_tolerance = 1e-6) {
  if (!inherits(prepared, "trsc_data") || !is.list(prepared) ||
      !identical(prepared$schema_version, 2L) ||
      !inherits(prepared$blocks, "trsc_blocks")) {
    stop("prepared must be a schema-version-2 trsc_data object.", call. = FALSE)
  }
  if (!is.numeric(kappa) || is.complex(kappa) || length(kappa) != 1L ||
      is.na(kappa) || !is.finite(kappa) || kappa < 0) {
    stop("kappa must be one finite nonnegative real number.", call. = FALSE)
  }
  if (!is.numeric(q) || is.complex(q) || length(q) != 1L ||
      is.na(q) || !is.finite(q) || !(q %in% c(1, 2))) {
    stop("q must be exactly 1 or 2.", call. = FALSE)
  }
  q <- as.integer(q)
  if (!is.logical(allow_inaccurate) || length(allow_inaccurate) != 1L ||
      is.na(allow_inaccurate)) {
    stop("allow_inaccurate must be one nonmissing logical value.", call. = FALSE)
  }
  if (!is.numeric(feasibility_tolerance) || is.complex(feasibility_tolerance) ||
      length(feasibility_tolerance) != 1L || is.na(feasibility_tolerance) ||
      !is.finite(feasibility_tolerance) || feasibility_tolerance < 0) {
    stop("feasibility_tolerance must be one finite nonnegative real number.", call. = FALSE)
  }
  solver <- match.arg(solver)

  blocks <- prepared$blocks
  all_periods <- blocks$periods
  if (is.null(periods)) {
    index <- seq_along(all_periods)
  } else {
    if (!is.numeric(periods) || is.complex(periods) || is.null(periods) ||
        length(periods) < 1L || anyNA(periods) || any(!is.finite(periods)) ||
        anyDuplicated(periods) || any(diff(periods) <= 0)) {
      stop("periods must be NULL or a strictly increasing finite numeric vector.", call. = FALSE)
    }
    index <- match(periods, all_periods)
    if (anyNA(index)) {
      stop("Every fitted period must belong to prepared$blocks$periods.", call. = FALSE)
    }
    if (length(index) > 1L && any(diff(index) != 1L)) {
      stop("Fitted periods must form one contiguous block.", call. = FALSE)
    }
  }
  if (length(index) < 3L || length(index) <= q) {
    stop("The fitting block must contain at least three dates and more dates than q.", call. = FALSE)
  }

  Y0 <- blocks$Y0[index, , drop = FALSE]
  Y1 <- as.numeric(blocks$Y1[index])
  fitted_periods <- all_periods[index]
  donor_names <- colnames(Y0)
  if (!is.matrix(Y0) || !is.numeric(Y0) || any(!is.finite(Y0)) ||
      any(!is.finite(Y1))) {
    stop("Prepared fitting outcomes must be finite numeric values.", call. = FALSE)
  }

  geometry <- .trsc_align_geometry(geometry, donor_names)
  problem <- .trsc_build_qp(Y1, Y0, geometry, kappa, q)

  result <- switch(
    solver,
    osqp = .trsc_solve_osqp(
      problem,
      solver_control = solver_control,
      allow_inaccurate = allow_inaccurate
    )
  )

  J <- ncol(Y0)
  n <- nrow(Y0)
  W <- matrix(result$x, nrow = J, ncol = n,
              dimnames = list(donor_names, as.character(fitted_periods)))

  equality_violation <- max(abs(colSums(W) - 1))
  lower_violation <- max(pmax(-W, 0))
  upper_violation <- max(pmax(W - 1, 0))
  max_violation <- max(equality_violation, lower_violation, upper_violation)
  if (!is.finite(max_violation) || max_violation > feasibility_tolerance) {
    stop(
      "Numerical solution exceeds feasibility_tolerance; maximum simplex ",
      "violation is ", format(max_violation, scientific = TRUE), ".",
      call. = FALSE
    )
  }

  fitted_values <- rowSums(Y0 * t(W))
  residuals <- Y1 - fitted_values
  roughness <- .trsc_temporal_energy(W, geometry, q)
  fit_loss <- sum(residuals^2)
  penalty <- kappa * roughness
  profile <- geometry$B %*% W
  colnames(profile) <- as.character(fitted_periods)

  structure(
    list(
      schema_version = 1L,
      call = match.call(),
      periods = fitted_periods,
      donor_names = donor_names,
      target_outcomes = Y1,
      donor_outcomes = Y0,
      geometry = geometry,
      specification = list(
        q = q,
        kappa = kappa,
        solver = solver,
        geometry = "supplied",
        reference_law_certified = FALSE
      ),
      predictive_path = profile,
      weights = W,
      fitted_values = fitted_values,
      residuals = residuals,
      objective = list(
        fit_loss = fit_loss,
        temporal_roughness = roughness,
        regularization_penalty = penalty,
        total = fit_loss + penalty
      ),
      solver = list(
        backend = "osqp",
        status = result$info$status,
        status_val = result$info$status_val,
        iterations = result$info$iter,
        primal_residual = result$info$prim_res,
        dual_residual = result$info$dual_res,
        objective_value = result$info$obj_val
      ),
      diagnostics = list(
        simplex_sum_violation = equality_violation,
        simplex_lower_violation = lower_violation,
        simplex_upper_violation = upper_violation,
        maximum_simplex_violation = max_violation,
        empirical_second_moment_rank = geometry$diagnostics$empirical_rank,
        predictive_dimension_empirical =
          geometry$diagnostics$predictive_dimension_empirical,
        simplex_null_dimension_empirical =
          geometry$diagnostics$simplex_null_dimension_empirical,
        representation_unique_under_empirical_geometry =
          geometry$diagnostics$representation_unique_empirical
      )
    ),
    class = c("trsc_fit", "trsc")
  )
}

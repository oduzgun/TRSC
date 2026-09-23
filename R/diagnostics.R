#' Extract empirical diagnostics from a TRSC fit
#'
#' @param object A `trsc_fit` object.
#'
#' @return A `trsc_diagnostics` object containing numerical fit, geometry,
#'   regularization and solver diagnostics.
#'
#' @details These quantities are descriptive diagnostics for the realized fit.
#'   They do not certify population coercivity, calibration complexity,
#'   continuation assumptions, or causal identification.
#'
#' @export
trsc_diagnostics <- function(object) {
  if (!inherits(object, "trsc_fit") || !is.list(object)) {
    stop("object must be a trsc_fit object.", call. = FALSE)
  }

  residuals <- as.numeric(object$residuals)
  structure(
    list(
      dates = length(object$periods),
      donors = length(object$donor_names),
      q = object$specification$q,
      kappa = object$specification$kappa,
      fit = list(
        sse = object$objective$fit_loss,
        rmse = sqrt(mean(residuals^2)),
        mae = mean(abs(residuals)),
        temporal_roughness = object$objective$temporal_roughness,
        regularization_penalty = object$objective$regularization_penalty,
        total_objective = object$objective$total
      ),
      geometry = list(
        empirical_second_moment_rank =
          object$diagnostics$empirical_second_moment_rank,
        predictive_dimension_empirical =
          object$diagnostics$predictive_dimension_empirical,
        simplex_null_dimension_empirical =
          object$diagnostics$simplex_null_dimension_empirical,
        representation_unique_under_empirical_geometry =
          object$diagnostics$representation_unique_under_empirical_geometry
      ),
      feasibility = list(
        simplex_sum_violation =
          object$diagnostics$simplex_sum_violation,
        simplex_lower_violation =
          object$diagnostics$simplex_lower_violation,
        simplex_upper_violation =
          object$diagnostics$simplex_upper_violation,
        maximum_simplex_violation =
          object$diagnostics$maximum_simplex_violation
      ),
      solver = object$solver,
      certification = list(
        population_coercivity = FALSE,
        calibration_complexity = FALSE,
        causal_identification = FALSE
      )
    ),
    class = "trsc_diagnostics"
  )
}

#' @export
print.trsc_diagnostics <- function(x, ...) {
  cat("TRSC empirical diagnostics\n")
  cat("  Dates: ", x$dates, "\n", sep = "")
  cat("  Donors: ", x$donors, "\n", sep = "")
  cat("  q: ", x$q, "\n", sep = "")
  cat("  Kappa: ", format(x$kappa), "\n", sep = "")
  cat("  RMSE: ", format(x$fit$rmse), "\n", sep = "")
  cat("  Temporal roughness: ", format(x$fit$temporal_roughness), "\n", sep = "")
  cat("  Empirical second-moment rank: ",
      x$geometry$empirical_second_moment_rank, "\n", sep = "")
  cat("  Empirical simplex null dimension: ",
      x$geometry$simplex_null_dimension_empirical, "\n", sep = "")
  cat("  Solver status: ", x$solver$status, "\n", sep = "")
  cat("  Maximum simplex violation: ",
      format(x$feasibility$maximum_simplex_violation, scientific = TRUE),
      "\n", sep = "")
  cat("  These diagnostics are not theorem or causal-identification certificates.\n")
  invisible(x)
}

#' Extract the fitted TRSC predictive path
#'
#' @param object A fitted TRSC object.
#' @param ... Additional arguments for methods.
#' @return Method-dependent predictive-path representation.
#' @importFrom stats fitted residuals weights
#' @export
predictive_path <- function(object, ...) {
  UseMethod("predictive_path")
}

#' @rdname predictive_path
#' @export
predictive_path.trsc_fit <- function(object, ...) {
  structure(
    list(
      periods = object$periods,
      donor_names = object$donor_names,
      representative = object$weights,
      metric_coordinates = object$metric_coordinates,
      geometry = object$geometry,
      diagnostics = list(
        empirical_second_moment_rank =
          object$diagnostics$empirical_second_moment_rank,
        simplex_null_dimension_empirical =
          object$diagnostics$simplex_null_dimension_empirical,
        representation_unique_under_empirical_geometry =
          object$diagnostics$representation_unique_under_empirical_geometry
      )
    ),
    class = "trsc_predictive_path"
  )
}

#' Print a TRSC predictive-path object
#'
#' @param x A `trsc_predictive_path` object.
#' @param ... Unused.
#' @export
print.trsc_predictive_path <- function(x, ...) {
  cat("TRSC fitted predictive relationship\n")
  cat("  Dates: ", length(x$periods), "\n", sep = "")
  cat("  Donors: ", length(x$donor_names), "\n", sep = "")
  cat("  Supplied-geometry rank: ",
      x$diagnostics$empirical_second_moment_rank, "\n", sep = "")
  cat("  Empirical representation null dimension: ",
      x$diagnostics$simplex_null_dimension_empirical, "\n", sep = "")
  cat("  Raw lift unique under supplied empirical geometry: ",
      if (isTRUE(x$diagnostics$representation_unique_under_empirical_geometry))
        "yes" else "no", "\n", sep = "")
  cat("  Note: metric coordinates are computational coordinates under the\n")
  cat("        supplied geometry, not a population-identification certificate.\n")
  invisible(x)
}

#' Extract one solver-selected donor representation
#'
#' @param object A `trsc_fit` object.
#' @param ... Unused.
#' @return A donor-by-time matrix. The matrix is a numerical simplex lift of
#'   the fitted predictive path. It is not automatically an identified
#'   structural donor-share path.
#' @export
weights.trsc_fit <- function(object, ...) {
  out <- object$weights
  attr(out, "representation") <- "solver-selected"
  attr(out, "unique_under_empirical_geometry") <-
    object$diagnostics$representation_unique_under_empirical_geometry
  out
}

#' @export
fitted.trsc_fit <- function(object, ...) {
  object$fitted_values
}

#' @export
residuals.trsc_fit <- function(object, ...) {
  object$residuals
}

#' Summarize a fitted TRSC model
#'
#' @param object A `trsc_fit` object.
#' @param ... Unused.
#' @return A `summary_trsc_fit` object.
#' @export
summary.trsc_fit <- function(object, ...) {
  structure(
    list(
      periods = range(object$periods),
      dates = length(object$periods),
      donors = length(object$donor_names),
      q = object$specification$q,
      kappa = object$specification$kappa,
      fit_loss = object$objective$fit_loss,
      temporal_roughness = object$objective$temporal_roughness,
      total_objective = object$objective$total,
      empirical_second_moment_rank =
        object$diagnostics$empirical_second_moment_rank,
      predictive_dimension_empirical =
        object$diagnostics$predictive_dimension_empirical,
      simplex_null_dimension_empirical =
        object$diagnostics$simplex_null_dimension_empirical,
      representation_unique_under_empirical_geometry =
        object$diagnostics$representation_unique_under_empirical_geometry,
      solver_status = object$solver$status,
      max_simplex_violation =
        object$diagnostics$maximum_simplex_violation
    ),
    class = "summary_trsc_fit"
  )
}

#' @export
print.summary_trsc_fit <- function(x, ...) {
  cat("TRSC pre-treatment fit\n")
  cat("  Fitted dates: ", x$dates, " (", x$periods[1], " to ",
      x$periods[2], ")\n", sep = "")
  cat("  Donors: ", x$donors, "\n", sep = "")
  cat("  Difference order q: ", x$q, "\n", sep = "")
  cat("  Kappa: ", format(x$kappa), "\n", sep = "")
  cat("  Empirical second-moment rank: ",
      x$empirical_second_moment_rank, "\n", sep = "")
  cat("  Predictive simplex dimension: ",
      x$predictive_dimension_empirical, "\n", sep = "")
  cat("  Simplex representation null dimension: ",
      x$simplex_null_dimension_empirical, "\n", sep = "")
  cat("  Representation unique under empirical geometry: ",
      if (isTRUE(x$representation_unique_under_empirical_geometry)) "yes" else "no",
      "\n", sep = "")
  cat("  Fit loss: ", format(x$fit_loss), "\n", sep = "")
  cat("  Temporal roughness: ", format(x$temporal_roughness), "\n", sep = "")
  cat("  Total objective: ", format(x$total_objective), "\n", sep = "")
  cat("  Solver status: ", x$solver_status, "\n", sep = "")
  cat("  Maximum simplex violation: ",
      format(x$max_simplex_violation, scientific = TRUE), "\n", sep = "")
  invisible(x)
}

#' @export
print.trsc_fit <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

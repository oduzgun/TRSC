#' Extract empirical diagnostics from a TRSC fit
#'
#' @param object A `trsc_fit`, identified/confidence set, or support result.
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
  if (inherits(object, "trsc_support")) {
    return(structure(list(backend=object$backend,status=object$status,
      primal_objective=object$primal_objective,certified_upper=object$upper,
      numerical_certificate=object$certificate,primal_feasibility=object$primal_feasibility,
      certified_gap=object$certified_gap,statistical_coverage=object$coverage,
      causal_assumptions_verified=FALSE),class="trsc_set_diagnostics"))
  }
  if (inherits(object, "trsc_identified_set")) {
    out<-list(type=object$type,H=object$H,structural_q=object$q,R=object$R,
              empirical_rank=object$geometry$diagnostics$empirical_rank,
              auxiliary_restrictions=length(object$auxiliary),
              causal_assumptions_verified=FALSE)
    if (inherits(object,"trsc_confidence_set")) {
      ac<-object$anchor_certificate;eps<-.trsc_Q(object$geometry_certificate$epsilon)
      t0<-2*.trsc_Q(ac$radius0)+2*eps
      b<-t0
      if(object$q==2) {
        tv<-2*.trsc_Q(ac$radius_velocity)+4*eps
        b<-.trsc_Q(.trsc_norm_up(c(t0+tv,if(object$H>=2)t0 else .trsc_Q(0))))
      }
      out$outer_radius_allowance<-.trsc_up(2*.trsc_Q(object$energy_allowance)+b)
      out$components<-ac$components;out$coverage<-object$coverage
      out$calibration_scope<-object$geometry_certificate$class_controlled
      out$anchor_propagation<-trsc_anchor_propagation(ac,seq_len(object$H))
      out$sandwich_scope<-"on the supplied common event; numerical endpoint enclosure error is separate"
    }
    return(structure(out,class="trsc_set_diagnostics"))
  }
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
print.trsc_set_diagnostics <- function(x, ...) {
  cat("TRSC set/optimization diagnostics\n")
  print(unclass(x))
  invisible(x)
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

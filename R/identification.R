# Profile/evaluation geometry. Coordinates are forward metric evaluations;
# anchors supplied as simplex representatives fix profiles, not future lifts.

.trsc_donor_vector <- function(x, geometry, label) {
  if (!is.numeric(x) || is.complex(x) || length(x) != geometry$donors || any(!is.finite(x)))
    stop(label, " must be a finite vector with one entry per donor.", call. = FALSE)
  if (!is.null(names(x))) {
    if (anyDuplicated(names(x)) || !setequal(names(x), geometry$donor_names))
      stop(label, " names must match donors exactly.", call. = FALSE)
    x <- x[geometry$donor_names]
  }
  as.numeric(x)
}

.trsc_simplex_rep <- function(x, geometry, label) {
  if (inherits(x, "trsc_readout")) x <- x$representative
  x <- .trsc_donor_vector(x, geometry, label)
  if (any(x < 0) || abs(sum(x) - 1) > 1e-7)
    stop(label, " must be a nonnegative simplex representative. Negative solver weights are not silently clipped.", call. = FALSE)
  z <- .trsc_Q(x)
  z / sum(z)
}

.trsc_deployment <- function(x, geometry) {
  if (is.vector(x)) x <- matrix(x, 1, dimnames = list(NULL, names(x)))
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
      ncol(x) != geometry$donors || nrow(x) < 1 || any(!is.finite(x)))
    stop("x must be a finite horizon-by-donor matrix.", call. = FALSE)
  if (!is.null(colnames(x))) {
    if (anyDuplicated(colnames(x)) || !setequal(colnames(x), geometry$donor_names))
      stop("x column names must match donors exactly.", call. = FALSE)
    x <- x[, geometry$donor_names, drop = FALSE]
  }
  colnames(x) <- geometry$donor_names
  x
}

.trsc_order <- function(q) {
  if (length(q) != 1 || !is.numeric(q) || is.na(q) || !q %in% c(1,2))
    stop("Structural q must be supplied explicitly as 1 or 2.", call. = FALSE)
  as.integer(q)
}

.trsc_anchors <- function(geometry, anchor, q, previous) {
  a0 <- .trsc_simplex_rep(anchor, geometry, "anchor")
  if (q == 1) return(list(a0))
  if (is.null(previous)) stop("q=2 requires the preceding anchor representative.", call. = FALSE)
  list(.trsc_simplex_rep(previous, geometry, "previous"), a0)
}

.trsc_continuation_map <- function(geometry, H, q, anchors) {
  J <- geometry$donors
  K <- .trsc_Q(kronecker(.trsc_difference_matrix(H + q, q), geometry$B))
  fixed <- seq_len(q * J)
  list(M = K[, -fixed, drop = FALSE],
       offset = .trsc_mm(K[, fixed, drop = FALSE], do.call(c, anchors)))
}

.trsc_dynamic_program <- function(object) {
  J <- object$geometry$donors; H <- object$H
  p <- .trsc_conic_problem(J * H, lapply(seq_len(H), function(h) (h-1L)*J + seq_len(J)))
  map <- .trsc_continuation_map(object$geometry, H, object$q, object$anchors)
  p <- .trsc_add_soc(p, map$M, map$offset, object$R)
  .trsc_add_auxiliary(p, object, future_start = 0L)
}

.trsc_add_auxiliary <- function(p, object, future_start) {
  J <- object$geometry$donors
  for (restriction in object$auxiliary) {
    for (h in seq_len(object$H)) {
      C <- restriction$maps[[h]]
      M <- .trsc_Q(matrix(0, nrow(C), p$n))
      M[, future_start + (h-1L)*J + seq_len(J)] <- .trsc_Q(C)
      if (restriction$error[h] == 0) {
        p <- .trsc_add_equality(p, M, restriction$values[[h]])
      } else {
        p <- .trsc_add_soc(p, M, -.trsc_Q(restriction$values[[h]]), restriction$error[h])
      }
    }
  }
  p
}

#' Evaluation fiber at a supplied predictive profile
#'
#' @param geometry A `trsc_geometry` object. Its forward factor defines the
#'   supplied prediction geometry; empirical geometry is not population truth.
#' @param x A finite named deployment donor vector.
#' @param profile Supplied forward coordinates, with length `nrow(geometry$B)`.
#' @param representative Alternative simplex lift of the profile. Supply
#'   exactly one of `profile` or `representative`.
#' @param control Named ECOS controls.
#' @return A `trsc_evaluation_set` with certified outer endpoint enclosures,
#'   raw optimizer values and compatible numerical lifts. Values are numerical
#'   evaluations of an identification problem, not sampling confidence limits.
#' @details The two LPs retain all lifts. Supplied representatives are converted
#'   to exact rational proportions (sum-one roundoff only, never negative-weight
#'   clipping). Their predictive profile is then held fixed exactly. An
#'   explicitly supplied coordinate vector defines its own exact constraints.
#' @export
trsc_evaluation_set <- function(geometry, x, profile = NULL, representative = NULL, control = list()) {
  .trsc_align_geometry(geometry, geometry$donor_names)
  x <- .trsc_donor_vector(x, geometry, "x")
  if (is.null(profile) == is.null(representative)) stop("Supply exactly one of profile or representative.", call. = FALSE)
  if (!is.null(representative)) profile <- .trsc_mm(geometry$B, .trsc_simplex_rep(representative, geometry, "representative"))
  if (length(profile) != nrow(geometry$B) || any(!is.finite(.trsc_num(profile)))) stop("profile dimension or values are invalid.", call. = FALSE)
  p <- .trsc_conic_problem(geometry$donors, list(seq_len(geometry$donors)))
  p <- .trsc_add_equality(p, geometry$B, profile)
  hi <- .trsc_solve_lp(p, x, control); lo <- .trsc_solve_lp(p, -x, control)
  empty <- hi$status == "structural_infeasible" || lo$status == "structural_infeasible"
  structure(list(interval = if (empty) c(lower=Inf, upper=-Inf) else c(lower=-lo$upper, upper=hi$upper),
                 status = if (empty) "structural_infeasible" else "computed",
                 upper_result = hi, lower_result = lo,
                 compatible_lifts = list(lower = lo$solution, upper = hi$solution),
                 geometry = geometry, x = x, profile = profile,
                 interpretation = "evaluation-fiber endpoint enclosure, not a confidence interval"),
            class = "trsc_evaluation_set")
}

#' Evaluation-fiber width and numerical point-identification diagnostic
#' @param object A `trsc_evaluation_set`.
#' @param tolerance Nonnegative numerical width tolerance. A TRUE result means
#'   the certified enclosure is no wider than this tolerance, not symbolic proof.
#' @return Width, or a logical diagnostic. Empty fibers return NA for the
#'   diagnostic and have status `structural_infeasible`.
#' @export
trsc_evaluation_width <- function(object) {
  if (!inherits(object, "trsc_evaluation_set")) stop("object must be an evaluation set.")
  if (object$status == "structural_infeasible") return(NA_real_)
  .trsc_up(.trsc_Q(object$interval[2]) - .trsc_Q(object$interval[1]))
}

#' @rdname trsc_evaluation_width
#' @export
trsc_point_identified <- function(object, tolerance = 1e-7) {
  tolerance <- .trsc_scalar(tolerance, "tolerance")
  trsc_evaluation_width(object) <= tolerance
}

#' Construct an assumption-indexed dynamic counterfactual set
#' @param geometry Supplied predictive geometry.
#' @param x Horizon-by-donor deployment matrix; row order is temporal order.
#' @param anchor Simplex representative (or frozen readout) of the cutoff profile.
#' @param q Explicit structural order, 1 or 2. Never inferred from a fitted q.
#' @param R Declared nonnegative root-energy sensitivity radius, not estimated.
#' @param previous For q=2, simplex representative of the preceding profile.
#' @param type `plugin` (default) or `population`; the latter is an explicit
#'   assertion that supplied geometry and causal anchors are population inputs.
#' @param control Named ECOS controls.
#' @return A lazy `trsc_identified_set` support oracle. This stores the entire
#'   feasible program; finitely many support calls do not enumerate that set.
#' @details Causal interpretation requires the untreated conditional-mean
#'   bridge and the declared continuation law. Fixing an anchor representative
#'   here fixes its profile only; future lifts can vary freely. No post-treatment
#'   treated outcome is used. The simplex ratio normalization is recorded.
#' @export
trsc_counterfactual_set <- function(geometry, x, anchor, q, R, previous = NULL,
                                    type = c("plugin", "population"), control = list()) {
  .trsc_align_geometry(geometry, geometry$donor_names)
  q <- .trsc_order(q); R <- .trsc_scalar(R, "R"); type <- match.arg(type)
  x <- .trsc_deployment(x, geometry)
  anchors <- .trsc_anchors(geometry, anchor, q, previous)
  structure(list(geometry = geometry, x = x, H = nrow(x), q = q, R = R,
                 anchors = anchors, type = type, auxiliary = list(), control = control,
                 interpretation = paste(type, "identified set under supplied causal bridge and continuation assumptions"),
                 anchor_normalization = "nonnegative input proportions normalized exactly; no clipping"),
            class = "trsc_identified_set")
}

#' Certified support enclosure of a dynamic set
#' @param object A `trsc_identified_set` or `trsc_confidence_set`.
#' @param a Finite horizon weights defining the linear functional.
#' @param control Optional replacement ECOS controls.
#' @return A `trsc_support` object: `upper` is a certified numerical upper
#'   enclosure; `primal_objective` is an approximate solver value, not an upper
#'   bound. Certified structural infeasibility is explicit with upper = -Inf.
#' @export
trsc_support <- function(object, a, control = object$control) {
  if (!inherits(object, "trsc_identified_set")) stop("object must be a TRSC set.", call. = FALSE)
  if (!is.numeric(a) || is.complex(a) || length(a) != object$H || any(!is.finite(a))) stop("a must have one finite weight per horizon.", call. = FALSE)
  confidence <- inherits(object, "trsc_confidence_set")
  p <- if (confidence) .trsc_confidence_program(object) else .trsc_dynamic_program(object)
  cfuture <- do.call(c, lapply(seq_len(object$H), function(h) .trsc_Q(a[h]) * .trsc_Q(object$x[h, ])))
  obj <- if (confidence) c(.trsc_Q(rep(0, object$q*object$geometry$donors)), cfuture) else cfuture
  ans <- .trsc_solve_socp(p, obj, control)
  ans$a <- as.numeric(a); ans$H <- object$H; ans$q <- object$q; ans$R <- object$R
  ans$uncertainty <- if (confidence) "outer confidence support" else paste(object$type, "identified support")
  ans$coverage <- if (confidence) object$coverage else NULL
  class(ans) <- c("trsc_support", "list")
  ans
}

#' Coordinate enclosures from a set support oracle
#' @param object A TRSC dynamic or confidence set.
#' @return A data frame with horizon, lower and upper enclosures. For a
#'   confidence set these are simultaneous projections on its common event.
#'   The rectangle generally loses the joint set's temporal restrictions.
#' @export
trsc_bounds <- function(object) {
  H <- object$H
  results <- lapply(seq_len(H), function(h) {
    a <- rep(0, H); a[h] <- 1
    list(lo = trsc_support(object,-a), hi = trsc_support(object,a))
  })
  out <- data.frame(horizon=seq_len(H),lower=vapply(results,function(x)-x$lo$upper,numeric(1)),
                    upper=vapply(results,function(x)x$hi$upper,numeric(1)))
  attr(out,"support_results") <- results
  attr(out,"interpretation") <- object$interpretation
  out
}

#' Minimum feasible continuation radius
#' @param geometry,anchor,q,previous,control As in [trsc_counterfactual_set()].
#' @param H Positive integer horizon.
#' @return `trsc_min_radius` with a certified lower bound, feasible upper
#'   bound and solver estimate of root energy. The upper witness is explicitly
#'   projected onto future simplices using exact rational proportions. Neither
#'   raw optimizer feasibility nor strict feasibility is assumed.
#' @export
trsc_min_radius <- function(geometry, anchor, q, H, previous = NULL, control = list()) {
  H <- .trsc_positive_integer(H,"H"); q <- .trsc_order(q)
  anchors <- .trsc_anchors(geometry, anchor, q, previous)
  if (q == 1L) return(structure(list(lower=0,upper=0,estimate=0,status="analytic_persistence"),class="trsc_min_radius"))
  J <- geometry$donors; n <- J*H
  map <- .trsc_continuation_map(geometry,H,q,anchors)
  maxnorm <- max(vapply(seq_len(J),function(j).trsc_norm_up(geometry$B[,j]),numeric(1)))
  cap <- .trsc_norm_up(rep(.trsc_Q(2^q)*.trsc_Q(maxnorm),H))
  p <- .trsc_conic_problem(n+1L,lapply(seq_len(H),function(h)(h-1L)*J+seq_len(J)),cap)
  G <- cbind(rbind(.trsc_Q(matrix(0,1,n)),-map$M), .trsc_Q(c(-1,rep(0,nrow(map$M)))))
  p$soc <- list(list(G=G,h=c(.trsc_Q(0),map$offset)))
  ans <- .trsc_solve_socp(p,c(rep(0,n),-1),control)
  w <- .trsc_Q(rep(1/J,n))
  if (length(ans$solution)==n+1 && all(is.finite(ans$solution))) {
    w <- .trsc_Q(pmax(ans$solution[seq_len(n)],0))
  }
  for (h in seq_len(H)) {
    ix <- (h-1L)*J+seq_len(J)
    if (sum(w[ix])==0) w[ix] <- .trsc_Q(rep(1,J))
    w[ix] <- w[ix]/sum(w[ix])
  }
  feasible_upper <- .trsc_norm_up(.trsc_mm(map$M,w)+map$offset)
  structure(list(lower=max(0,-ans$upper),upper=feasible_upper,estimate=-ans$primal_objective,
                 status=ans$status,solver=ans,feasible_witness=w),class="trsc_min_radius")
}

#' Sensitivity curve over declared structural radii
#' @param object A dynamic or confidence set.
#' @param radii A finite nonnegative radius vector; no radius is estimated.
#' @param a Horizon functional weights.
#' @return Data frame of certified endpoint enclosures, retaining support results.
#' @export
trsc_sensitivity <- function(object, radii, a) {
  if (!is.numeric(radii) || !length(radii) || any(!is.finite(radii)) || any(radii<0)) stop("radii must be finite and nonnegative.")
  results <- lapply(radii,function(R) {
    z<-object;z$R<-R
    list(lo=trsc_support(z,-a),hi=trsc_support(z,a))
  })
  out<-data.frame(R=radii,lower=vapply(results,function(x)-x$lo$upper,numeric(1)),
                  upper=vapply(results,function(x)x$hi$upper,numeric(1)))
  attr(out,"support_results")<-results
  out
}

#' Refine a dynamic set using the same latent witnesses
#' @param object A dynamic set.
#' @param maps List of horizon-specific measurement matrices, columns in donor
#'   order. These must be scientifically justified augmented response maps,
#'   for example R_h P_D, not arbitrary restrictions on raw donor shares.
#' @param values List of corresponding untreated measurement vectors.
#' @param source Nonempty explanation of shared composition, identified response
#'   maps and treatment exclusion. Merely observing an auxiliary outcome is not
#'   sufficient to establish these restrictions.
#' @param error Nonnegative deterministic norm-error bounds, one per horizon
#'   (scalar recycled); zero imposes exact equations. A noisy realization cannot
#'   be used as an exact mean equation without an appropriate model.
#' @return Refined set. Measurement restrictions intersect the latent feasible
#'   path before evaluation/projection, retaining donor identity.
#' @export
trsc_auxiliary_refine <- function(object,maps,values,source,error=0) {
  if (!inherits(object,"trsc_identified_set")) stop("object must be a TRSC set.")
  if (!is.character(source)||length(source)!=1||is.na(source)||!nzchar(source)) stop("source must describe the structural measurement restriction.")
  if (!is.list(maps)||!is.list(values)||length(maps)!=object$H||length(values)!=object$H) stop("One map and value vector are required per horizon.")
  error<-rep_len(error,object$H)
  if (!is.numeric(error)||any(!is.finite(error))||any(error<0)) stop("error must be nonnegative.")
  for(h in seq_len(object$H)) {
    M<-maps[[h]];v<-values[[h]]
    if(!is.matrix(M)||!is.numeric(M)||is.complex(M)||ncol(M)!=object$geometry$donors||any(!is.finite(M))||
       !is.numeric(v)||is.complex(v)||length(v)!=nrow(M)||any(!is.finite(v))) stop("Auxiliary dimensions or values are invalid.")
    if(!is.null(colnames(M))) {
      if(anyDuplicated(colnames(M))||!setequal(colnames(M),object$geometry$donor_names)) stop("Auxiliary donor names differ.")
      maps[[h]]<-M[,object$geometry$donor_names,drop=FALSE]
    }
  }
  object$auxiliary<-c(object$auxiliary,list(list(maps=maps,values=values,error=error,source=source)))
  object
}

#' @export
print.trsc_identified_set <- function(x,...) {
  cat("TRSC ",x$type," set: H=",x$H,", structural q=",x$q,", R=",x$R,"\n",sep="")
  cat("  ",x$interpretation,"\n",sep="")
  cat("  Future lifts remain free subject to profile and supplied measurement constraints.\n")
  invisible(x)
}

#' @export
print.trsc_support <- function(x,...) {
  cat("TRSC ",x$uncertainty,"\n  Certified upper enclosure: ",format(x$upper),
      "\n  Solver status: ",x$status,"\n  Certificate: ",x$certificate$method,"\n",sep="")
  invisible(x)
}

#' @export
print.trsc_evaluation_set <- function(x,...) {
  cat("TRSC evaluation fiber: ",x$status,"\n",sep="");print(x$interval)
  cat("  Numerical endpoint enclosure, not sampling inference.\n");invisible(x)
}

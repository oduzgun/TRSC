.trsc_certificate_metadata <- function(failure_probability,source,provenance,event_id,conditioning) {
  failure_probability <- .trsc_scalar(failure_probability,"failure_probability",0,1)
  if(failure_probability>=1) stop("failure_probability must be less than one.")
  for(v in list(source,event_id)) if(!is.character(v)||length(v)!=1||is.na(v)||!nzchar(v)) stop("source and event_id must be nonempty strings.")
  provenance<-match.arg(provenance,c("maintained_bound","analytically_known","independently_certified","diagnostic_only"))
  conditioning<-match.arg(conditioning,c("fixed_or_unconditional","conditional_on_deployment"))
  list(failure_probability=failure_probability,source=source,provenance=provenance,
       event_id=event_id,conditioning=conditioning)
}

#' Explicit causal-anchor uncertainty certificate
#' @param geometry Empirical geometry with the intended donor order.
#' @param center Nonnegative simplex representative (or readout) of the fitted
#'   level center. Sum-one roundoff is normalized as exact rational proportions.
#'   Supplied bounds must apply to this recorded center, including any numerical
#'   correction. Negative solver weights are rejected, never silently clipped.
#' @param q Boundary order, 1 or 2. A q2 certificate may also support q1 inference.
#' @param components Named nonnegative vector with `estimation`, `fitting_gap`,
#'   `causal_gap`, `numerical`. For q1 these bound prediction-norm errors; for
#'   q2 they bound the joint W norm of level and velocity errors.
#' @param failure_probability Probability budget for the common path/anchor
#'   event including any fitting exceptional event. It is not charged again for
#'   each readout or each horizon.
#' @param source Explanation or reference supporting the numerical bounds.
#' @param previous Preceding readout representative for q2; velocity center is
#'   center minus previous in predictive geometry.
#' @param W Positive definite symmetric 2-by-2 joint weighting matrix for q2.
#' @param provenance One of `maintained_bound`, `analytically_known`,
#'   `independently_certified`, `diagnostic_only`.
#' @param event_id Identifier for a genuinely common event. Matching identifiers
#'   deduplicate probability only with matching probability and declared scope.
#' @param candidate_id Identifier of this frozen candidate.
#' @param covered_candidates Candidate identifiers jointly covered by this event.
#' @param conditioning `fixed_or_unconditional`, or
#'   `conditional_on_deployment` only if the certificate actually holds conditionally.
#' @return A `trsc_anchor_certificate`. Construction validates the input contract,
#'   not the underlying economic/statistical assumptions. Unknown risk envelopes,
#'   coercivity or smoothness constants are never estimated implicitly.
#' @details Regularization alone does not bound the fitting-to-cutoff gap.
#'   Predictive error alone does not establish the causal anchor bridge. For q2
#'   both marginal radii inherit from the single joint event, with exact rational
#'   inversion of W (never inversion of predictive geometry).
#' @export
trsc_anchor_certificate <- function(geometry,center,q,components,failure_probability,source,
                                    previous=NULL,W=diag(2),provenance="maintained_bound",
                                    event_id="anchor",candidate_id="candidate",
                                    covered_candidates=candidate_id,
                                    conditioning="fixed_or_unconditional") {
  q<-.trsc_order(q);.trsc_align_geometry(geometry,geometry$donor_names)
  required<-c("estimation","fitting_gap","causal_gap","numerical")
  if(!is.numeric(components)||is.complex(components)||any(!is.finite(components))||any(components<0)||
     length(components)!=4||anyDuplicated(names(components))||!setequal(names(components),required))
    stop("components must name estimation, fitting_gap, causal_gap and numerical bounds.")
  components<-components[required]
  if(!is.character(candidate_id)||length(candidate_id)!=1||is.na(candidate_id)||!nzchar(candidate_id)||
     !is.character(covered_candidates)||anyNA(covered_candidates)||anyDuplicated(covered_candidates)||
     !candidate_id %in% covered_candidates) stop("Invalid candidate/event scope.")
  meta<-.trsc_certificate_metadata(failure_probability,source,provenance,event_id,conditioning)
  anchors<-.trsc_anchors(geometry,center,q,previous)
  rho<-.trsc_up(sum(.trsc_Q(components)))
  out<-c(meta,list(geometry=geometry,q=q,anchors=anchors,components=components,
                    candidate_id=candidate_id,covered_candidates=covered_candidates,
                    radius0=rho,radius_velocity=NULL,joint_radius=NULL,W=NULL))
  if(q==2) {
    if(!is.matrix(W)||!is.numeric(W)||!identical(dim(W),c(2L,2L))||any(!is.finite(W))||
       W[1,2]!=W[2,1]) stop("W must be a finite symmetric 2-by-2 matrix.")
    Wq<-.trsc_Q(W);detW<-Wq[1,1]*Wq[2,2]-Wq[1,2]^2
    if(Wq[1,1]<=0||detW<=0) stop("W must be positive definite.")
    inv<-rbind(c(Wq[2,2],-Wq[1,2]),c(-Wq[1,2],Wq[1,1]))/c(detW)
    # Upper square roots verified by exact squared comparisons.
    sqrt_up<-function(v) .trsc_sqrt_up(v)
    out$radius0<-sqrt_up(.trsc_Q(rho)^2*inv[1,1])
    out$radius_velocity<-sqrt_up(.trsc_Q(rho)^2*inv[2,2])
    out$joint_radius<-rho;out$W<-W;out$W_inverse_exact<-inv
    U<-chol(W);out$W_root<-U
    residual<-.trsc_mm(t(U),U)-Wq
    delta<-max(sum(abs(residual[1,])),sum(abs(residual[2,])))
    out$joint_computational_radius<-sqrt_up(.trsc_Q(rho)^2+delta*(.trsc_Q(out$radius0)^2+.trsc_Q(out$radius_velocity)^2))
    out$W_roundoff_allowance<-out$joint_computational_radius-rho
    out$calibration_constant<-.trsc_norm_up(abs(.trsc_Q(U[,1]))+2*abs(.trsc_Q(U[,2])))
  }
  structure(out,class="trsc_anchor_certificate")
}

.trsc_sqrt_up <- function(q) {
  q<-.trsc_Q(q)
  if(q<0) stop("Negative square-root argument.")
  v<-sqrt(as.numeric(.trsc_num(q)))
  if(!is.finite(v)) stop("Square root exceeds numeric range.")
  while(.trsc_Q(v)^2<q) v<-v+max(abs(v)*4*.Machine$double.eps,.Machine$double.xmin)
  v
}

#' Certificate for future-path chord calibration
#' @param geometry The actual learned geometry used by the continuation program.
#' @param epsilon Nonnegative certified absolute norm deviation on all simplex
#'   chords. Fitting-local calibration or an empirical rank is insufficient.
#' @param failure_probability,source,provenance,event_id,conditioning As in
#'   [trsc_anchor_certificate()].
#' @return A `trsc_geometry_certificate`; no population quantity is inferred from
#'   the rank of the calibration matrix. This is an explicit supplied certificate.
#' @export
trsc_geometry_certificate <- function(geometry,epsilon,failure_probability,source,
                                      provenance="maintained_bound",event_id="future_geometry",
                                      conditioning="fixed_or_unconditional") {
  .trsc_align_geometry(geometry,geometry$donor_names)
  meta<-.trsc_certificate_metadata(failure_probability,source,provenance,event_id,conditioning)
  structure(c(meta,list(geometry=geometry,epsilon=.trsc_scalar(epsilon,"epsilon"),
                        class_controlled="all simplex chords (future paths), not fitting-local errors")),
            class="trsc_geometry_certificate")
}

#' Boundary error propagation by horizon
#' @param certificate An anchor certificate.
#' @param horizons Nonnegative integer horizons.
#' @return Data frame of level-plus-velocity and, for q2, joint-ellipse upper
#'   bounds. These are profile errors; deployment amplification is separate.
#' @export
trsc_anchor_propagation <- function(certificate,horizons) {
  if(!inherits(certificate,"trsc_anchor_certificate")) stop("An anchor certificate is required.")
  if(!is.numeric(horizons)||any(!is.finite(horizons))||any(horizons<0)||any(horizons!=floor(horizons))) stop("Invalid horizons.")
  a<-certificate
  marginal<-vapply(horizons,function(h).trsc_up(.trsc_Q(a$radius0)+if(a$q==2).trsc_Q(h)*.trsc_Q(a$radius_velocity) else .trsc_Q(0)),numeric(1))
  joint<-if(a$q==2) vapply(horizons,function(h) {
    v<-.trsc_Q(c(1,h));.trsc_sqrt_up(.trsc_Q(a$joint_radius)^2*sum(v*.trsc_mm(a$W_inverse_exact,v)))
  },numeric(1)) else marginal
  data.frame(horizon=horizons,level_plus_velocity=marginal,joint=joint,bound=pmin(marginal,joint))
}

#' Form a common-event anchor certificate from a supplied risk envelope
#' @param geometry,center,q,source,previous,W,provenance,event_id,candidate_id,
#'   covered_candidates,conditioning As in [trsc_anchor_certificate()].
#' @param path_risk_bound Supplied bound Bbar for E[1_good ||error||^2]/n,
#'   uniform over admissible approximate minimizers. Not an oracle point estimate.
#' @param alpha_path Markov failure allocation on the common path event.
#' @param fitting_failure Additional probability that the fitting good event fails.
#' @param readout_weights One (q1) or two (q2) rows of frozen convex readout
#'   weights over the fitting dates: level, then preceding level.
#' @param fitting_gap,causal_gap,numerical Maintained bounds in the scalar or
#'   joint W norm. Numerical allowance must include any difference between the
#'   supplied normalized centers and the theoretical readouts.
#' @return Anchor certificate retaining the common path-event provenance.
#' @export
trsc_anchor_from_risk <- function(geometry,center,q,path_risk_bound,alpha_path,
                                  fitting_failure,readout_weights,fitting_gap,causal_gap,
                                  numerical,source,previous=NULL,W=diag(2),
                                  provenance="maintained_bound",event_id="common_path",
                                  candidate_id="candidate",covered_candidates=candidate_id,
                                  conditioning="fixed_or_unconditional") {
  q<-.trsc_order(q);B<-.trsc_scalar(path_risk_bound,"path_risk_bound")
  alpha_path<-.trsc_scalar(alpha_path,"alpha_path",.Machine$double.xmin,1)
  fitting_failure<-.trsc_scalar(fitting_failure,"fitting_failure",0,1)
  L<-readout_weights
  if(is.vector(L)) L<-matrix(L,1)
  if(!is.matrix(L)||!is.numeric(L)||nrow(L)!=q||any(!is.finite(L))||any(L<0)||any(abs(rowSums(L)-1)>1e-10)) stop("readout_weights must contain q convex rows.")
  n<-ncol(L);Lq<-.trsc_Q(L)
  if(q==2) Lq<-rbind(Lq[1,,drop=FALSE],Lq[1,,drop=FALSE]-Lq[2,,drop=FALSE])
  if(q==1) op2<-sum(Lq^2) else {
    # Largest eigenvalue of W^(1/2) L L' W^(1/2), via the similar
    # 2-by-2 matrix W L L'. Trace/determinant are exact rational quantities.
    M<-.trsc_mm(.trsc_Q(W),.trsc_mm(Lq,t(Lq)))
    tr<-M[1,1]+M[2,2];det<-M[1,1]*M[2,2]-M[1,2]*M[2,1]
    op2<-(tr+.trsc_Q(.trsc_sqrt_up(max(.trsc_Q(0),tr^2-4*det))))/2
  }
  estimation<-.trsc_sqrt_up(.trsc_Q(n)*.trsc_Q(B)*op2/.trsc_Q(alpha_path))
  out<-trsc_anchor_certificate(geometry,center,q,
    c(estimation=estimation,fitting_gap=fitting_gap,causal_gap=causal_gap,numerical=numerical),
    .trsc_up(.trsc_Q(alpha_path)+.trsc_Q(fitting_failure)),source,previous,W,provenance,
    event_id,candidate_id,covered_candidates,conditioning)
  out$common_path_event<-list(n=n,path_risk_bound=B,alpha_path=alpha_path,
     fitting_failure=fitting_failure,squared_path_radius=.trsc_up(.trsc_Q(n)*.trsc_Q(B)/.trsc_Q(alpha_path)),
     uniform_over_approximate_minimizers="required by supplied envelope",
     level_and_velocity_share_event=TRUE)
  out
}

#' @export
print.trsc_anchor_certificate <- function(x,...) {
  cat("TRSC anchor certificate (",x$provenance,"): q=",x$q,"\n",sep="")
  print(x$components);cat("  Common-event failure budget: ",x$failure_probability,"\n",sep="")
  cat("  Supplied assumptions are not statistically verified by this constructor.\n");invisible(x)
}

#' @export
print.trsc_geometry_certificate <- function(x,...) {
  cat("TRSC future-chord certificate: epsilon=",x$epsilon,"; budget=",x$failure_probability,"\n",sep="")
  cat("  Provenance: ",x$provenance,"\n",sep="");invisible(x)
}

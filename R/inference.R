.trsc_coverage <- function(ac,gc,library,mixing_error,conditional) {
  if(is.null(library)) library<-list(ac)
  if(!is.list(library)||!length(library)||!all(vapply(library,inherits,logical(1),"trsc_anchor_certificate"))) stop("library must contain anchor certificates.")
  ids<-vapply(library,`[[`,character(1),"candidate_id")
  if(anyDuplicated(ids)||!ac$candidate_id %in% ids||!identical(ac,library[[match(ac$candidate_id,ids)]])) stop("Selected anchor must match its frozen library certificate.")
  for(c in library) {
    sharing<-ids[vapply(library,function(z)identical(z$event_id,c$event_id),logical(1))]
    if(!all(sharing %in% c$covered_candidates)) stop("Shared event does not declare all covered candidates.")
  }
  certs<-c(library,list(gc))
  if(any(vapply(certs,function(z)z$provenance=="diagnostic_only",logical(1)))) stop("A diagnostic is not an inference certificate.")
  if(conditional&&!all(vapply(certs,function(z)z$conditioning=="conditional_on_deployment",logical(1)))) stop("Conditional coverage requires conditional certificates for every input.")
  if(gc$event_id %in% vapply(library,`[[`,character(1),"event_id")) stop("Geometry and anchor events need distinct identifiers.")
  events<-split(certs,vapply(certs,`[[`,character(1),"event_id"))
  probs<-vapply(events,function(z) {
    p<-vapply(z,`[[`,numeric(1),"failure_probability")
    if(length(unique(p))!=1) stop("A shared event has inconsistent failure probabilities.")
    p[1]
  },numeric(1))
  total<-.trsc_up(sum(.trsc_Q(probs))+.trsc_Q(mixing_error))
  if(total>=1) stop("Total failure allocation must be less than one.")
  list(failure_probability=total,nominal_coverage=1-total,events=probs,
       mixing_error=mixing_error,candidates=ids,selected=ac$candidate_id,
       simultaneous_horizons=TRUE,simultaneous_radii=TRUE,simultaneous_support_directions=TRUE,
       structural_order_validated=FALSE,conditional_on_deployment=conditional,
       scope=if(conditional) "conditional under supplied conditional certificates" else "fixed deployment experiment or unconditional random deployment; not conditional coverage")
}

#' Honest outer confidence-set support oracle
#' @param x Named horizon-by-donor deployment matrix.
#' @param anchor_certificate Explicit causal-boundary certificate.
#' @param geometry_certificate Explicit future-chord certificate.
#' @param q Explicit structural continuation order, never inherited from fitting.
#' @param R Nonnegative structural root-energy radius.
#' @param library Optional frozen list of candidate anchor certificates. The
#'   selected certificate must be in this list. Distinct event probabilities
#'   are added without independence assumptions; a common path event is counted
#'   once. Without a library, inference is for the supplied fixed candidate.
#' @param mixing_error Additional supplied coupling-error bound.
#' @param conditional Request conditional-on-deployment coverage; rejected unless
#'   all supplied certificates explicitly have that scope.
#' @param control ECOS numerical controls.
#' @return A `trsc_confidence_set` lazy joint outer set. Coverage is conditional
#'   on the validity of the numerical/statistical certificates and causal model.
#' @details Uses forward seminorms only. Singular learned geometry is allowed.
#'   Simultaneous radii and directions follow from one deterministic inclusion;
#'   finite candidates require probability allocation. Equal Markov allocations
#'   can cost sqrt(M/alpha), not merely a logarithm. No post-selection refit is
#'   covered. Predictive selection cannot validate structural q or R.
#' @export
trsc_confidence_set <- function(x,anchor_certificate,geometry_certificate,q,R,
                                library=NULL,mixing_error=0,conditional=FALSE,control=list()) {
  ac<-anchor_certificate;gc<-geometry_certificate
  if(!inherits(ac,"trsc_anchor_certificate")) stop("anchor_certificate is required.",call.=FALSE)
  if(!inherits(gc,"trsc_geometry_certificate")) stop("geometry_certificate is required.",call.=FALSE)
  if(!identical(ac$geometry$B,gc$geometry$B)||!identical(ac$geometry$donor_names,gc$geometry$donor_names)) stop("Anchor and geometry certificates refer to different geometry.")
  q<-.trsc_order(q);R<-.trsc_scalar(R,"R")
  if(q>ac$q) stop("q2 inference requires a joint level-velocity certificate.")
  if(!is.logical(conditional)||length(conditional)!=1||is.na(conditional)) stop("conditional must be TRUE or FALSE.")
  mixing_error<-.trsc_scalar(mixing_error,"mixing_error",0,1)
  x<-.trsc_deployment(x,gc$geometry);H<-nrow(x)
  coverage<-.trsc_coverage(ac,gc,library,mixing_error,conditional)
  allowance<-.trsc_norm_up(rep(.trsc_Q(2^(q-1))*.trsc_Q(gc$epsilon),H))
  structure(list(geometry=gc$geometry,x=x,H=H,q=q,R=R,type="confidence",
    anchor_certificate=ac,geometry_certificate=gc,coverage=coverage,
    energy_allowance=allowance,auxiliary=list(),control=control,
    interpretation="honest outer set under supplied causal, transfer and statistical certificates"),
    class=c("trsc_confidence_set","trsc_identified_set"))
}

.trsc_confidence_program <- function(object) {
  J<-object$geometry$donors;H<-object$H;q<-object$q;F<-.trsc_Q(object$geometry$B)
  ac<-object$anchor_certificate;eps<-.trsc_Q(object$geometry_certificate$epsilon)
  n<-J*(H+q)
  p<-.trsc_conic_problem(n,lapply(seq_len(H+q),function(h)(h-1L)*J+seq_len(J)))
  E0<-.trsc_Q(matrix(0,J,n));E0[,((q-1L)*J+seq_len(J))]<-.trsc_Q(diag(J))
  w0<-ac$anchors[[length(ac$anchors)]]
  p<-.trsc_add_soc(p,.trsc_mm(F,E0),-.trsc_mm(F,w0),.trsc_Q(ac$radius0)+eps)
  if(q==2) {
    Ev<-E0;Ev[,seq_len(J)]<- -.trsc_Q(diag(J))
    wv<-w0-ac$anchors[[1]]
    p<-.trsc_add_soc(p,.trsc_mm(F,Ev),-.trsc_mm(F,wv),.trsc_Q(ac$radius_velocity)+2*eps)
    U<-.trsc_Q(ac$W_root)
    M<-do.call(rbind,lapply(1:2,function(i).trsc_mm(F,c(U[i,1])*E0+c(U[i,2])*Ev)))
    off<-do.call(c,lapply(1:2,function(i)-.trsc_mm(F,c(U[i,1])*w0+c(U[i,2])*wv)))
    p<-.trsc_add_soc(p,M,off,.trsc_Q(ac$joint_computational_radius)+.trsc_Q(ac$calibration_constant)*eps)
  }
  K<-.trsc_Q(kronecker(.trsc_difference_matrix(H+q,q),object$geometry$B))
  p<-.trsc_add_soc(p,K,rep(0,nrow(K)),.trsc_Q(object$R)+.trsc_Q(object$energy_allowance))
  .trsc_add_auxiliary(p,object,future_start=q*J)
}

#' Certified upper endpoint of an outer confidence set
#' @param object A `trsc_confidence_set`.
#' @param a Horizon functional weights.
#' @param control ECOS controls.
#' @return `trsc_support`, with the certified endpoint separate from the raw
#'   primal optimizer. Lower endpoint is minus the upper endpoint at -a.
#' @export
trsc_confidence_support <- function(object,a,control=object$control) {
  if(!inherits(object,"trsc_confidence_set")) stop("A confidence set is required.")
  trsc_support(object,a,control)
}

#' Conditional-mean treatment-effect interval
#' @param object A TRSC outer confidence set.
#' @param a Horizon functional weights.
#' @param treated_mean Known treated conditional-mean functional a'mI, not an
#'   observed treated realization. Supply exactly one treated-mean input.
#' @param treated_mean_interval Certified interval for a'mI.
#' @param mean_failure_probability Failure allocation for that interval; required
#'   when it is supplied. Added by Bonferroni without independence assumptions.
#' @param mean_type Must equal `conditional`; realized individual effects and
#'   generic ATT estimands are deliberately not inferred by this interface.
#' @return A `trsc_effect_interval` for the declared conditional-mean effect.
#' @export
trsc_effect_interval <- function(object,a,treated_mean=NULL,treated_mean_interval=NULL,
                                  mean_failure_probability=NULL,mean_type="conditional") {
  if(!inherits(object,"trsc_confidence_set")) stop("object must be a confidence set.")
  if(!identical(mean_type,"conditional")) stop("The estimand requires a treated conditional mean, not a realized outcome.")
  if(is.null(treated_mean)==is.null(treated_mean_interval)) stop("Supply exactly one treated conditional-mean input.")
  if(!is.null(treated_mean)) {
    .trsc_scalar(treated_mean,"treated_mean",-Inf,Inf);mi<-rep(treated_mean,2);alphaI<-0
  } else {
    mi<-treated_mean_interval
    if(!is.numeric(mi)||length(mi)!=2||any(!is.finite(mi))||mi[1]>mi[2]) stop("Invalid treated-mean interval.")
    alphaI<-.trsc_scalar(mean_failure_probability,"mean_failure_probability",0,1)
  }
  upper<-trsc_confidence_support(object,a);negative<-trsc_confidence_support(object,-a)
  if(!is.finite(upper$upper)||!is.finite(negative$upper)) stop("Counterfactual model is infeasible; no effect interval is reported.")
  failure<-.trsc_up(.trsc_Q(object$coverage$failure_probability)+.trsc_Q(alphaI))
  if(failure>=1) stop("Combined failure allocation must be less than one.")
  interval<-c(lower=-.trsc_up(.trsc_Q(upper$upper)-.trsc_Q(mi[1])),
              upper=.trsc_up(.trsc_Q(mi[2])+.trsc_Q(negative$upper)))
  structure(list(interval=interval,failure_probability=failure,a=a,
                 estimand="a'(treated conditional mean - untreated conditional mean)",
                 support=list(upper=upper,negative=negative),realized_individual_effect=FALSE),
            class="trsc_effect_interval")
}

#' Deployment modulus on the full predictive hull
#' @param geometry Supplied predictive geometry.
#' @param x Named deployment vector.
#' @param epsilon Nonnegative predictive separation radius.
#' @param control ECOS controls.
#' @return Certified upper enclosure of the maximum evaluation difference
#'   between two simplex lifts at predictive distance at most epsilon. At zero
#'   this detects evaluation ambiguity. It is not by itself a universal minimax
#'   bound: a lower-bound application needs an admissible two-point experiment.
#' @export
trsc_deployment_modulus <- function(geometry,x,epsilon,control=list()) {
  epsilon<-.trsc_scalar(epsilon,"epsilon");x<-.trsc_donor_vector(x,geometry,"x");J<-geometry$donors
  p<-.trsc_conic_problem(2*J,list(seq_len(J),J+seq_len(J)))
  p<-.trsc_add_soc(p,cbind(.trsc_Q(geometry$B),-.trsc_Q(geometry$B)),rep(0,nrow(geometry$B)),epsilon)
  out<-.trsc_solve_socp(p,c(x,-x),control)
  out$epsilon<-epsilon;out$uncertainty<-"deployment modulus enclosure"
  class(out)<-c("trsc_support","list");out
}

#' @export
print.trsc_confidence_set <- function(x,...) {
  cat("TRSC honest outer set (under supplied certificates)\n")
  cat("  H=",x$H,", structural q=",x$q,", R=",x$R,"\n",sep="")
  cat("  Failure allocation: ",x$coverage$failure_probability,"\n",sep="")
  cat("  ",x$coverage$scope,"\n",sep="");invisible(x)
}

#' @export
print.trsc_effect_interval <- function(x,...) {
  cat("TRSC conditional-mean effect interval\n");print(x$interval)
  cat("  Failure allocation: ",x$failure_probability,"\n",sep="");invisible(x)
}

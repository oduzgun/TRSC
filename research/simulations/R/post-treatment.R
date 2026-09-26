# Research-only post-treatment harness. No exported package function receives
# this generator's truth object. Source after library(TRSC).

trsc_sim_post_ids <- function() c("persistence", "affine", "hull_exit",
  "exact_ambiguity", "near_redundancy", "weak_excitation", "full_rank",
  "curved", "abrupt")

trsc_sim_post_scenario <- function(scenario = "persistence", seed = 101L,
    T = 120L, J = 3L, H = 5L, rank = 1L, noise_sd = .2, gap = 10L,
    velocity = .05/T, amplification = 1, near_distance = 1e-5,
    dependence = 0) {
  if (!scenario %in% trsc_sim_post_ids()) stop("Unknown post-treatment scenario.")
  if (T < 30 || J < 2 || H < 1 || rank < 1 || rank > J || gap < 3 ||
      gap >= T/2 || noise_sd < 0 || dependence < 0 || dependence >= 1 ||
      amplification < 0 || near_distance < 0 || near_distance > 1)
    stop("Invalid post-treatment design cell.")
  set.seed(seed)
  donors <- paste0("D", seq_len(J))
  r <- if (scenario == "full_rank") J else rank
  F <- matrix(0, r, J, dimnames = list(NULL, donors))
  F[1, ] <- c(0, 1, if(J > 2) seq(.2, .8, length.out = J-2))
  if (r > 1) for (j in 2:r) F[j, j] <- 1
  if (r == J) F <- diag(J)
  colnames(F) <- donors
  F <- sweep(F, 2, pmax(1, sqrt(colSums(F^2))), "/")
  if (scenario == "exact_ambiguity") F[, 2] <- F[, 1]
  if (scenario %in% c("near_redundancy", "weak_excitation")) {
    F[, 2] <- F[, 1]; F[1, 2] <- F[1, 1] + near_distance
  }
  # Known population second moment of the latent sign state is I.
  population_geometry <- trsc_geometry(sqrt(r) * F)
  F <- population_geometry$B
  state <- matrix(sample(c(-1, 1), T*r, replace = TRUE), T, r)
  if (dependence > 0) for (t in 2:T) {
    if (runif(1) < dependence) state[t, ] <- state[t-1, ]
  }
  X <- state %*% F; colnames(X) <- donors
  b <- rep(.4, T + H)
  if (scenario == "affine") b <- .4 + ((1:(T+H))-T)*velocity
  if (scenario == "hull_exit") { b[seq_len(T-1)] <- 0; b[T:(T+H)] <- 1 }
  if (scenario == "curved") b <- .5 + .2*sin((1:(T+H))/T*2*pi)
  if (scenario == "abrupt") b[(T+1):(T+H)] <- .9
  if (any(b < 0 | b > 1)) stop("Requested relationship path leaves the simplex.")
  W <- matrix(0, J, T+H, dimnames = list(donors, NULL))
  W[1, ] <- 1-b; W[2, ] <- b
  mu_pre <- rowSums(X * t(W[, 1:T, drop = FALSE]))
  error <- rnorm(T, sd = noise_sd)
  if (dependence > 0) for (t in 2:T)
    error[t] <- dependence*error[t-1] + sqrt(1-dependence^2)*error[t]
  panel <- expand.grid(unit = c("Target", donors), period = seq_len(T),
                       stringsAsFactors = FALSE)
  panel$outcome <- c(rbind(mu_pre+error, t(X)))
  deployment <- matrix(rep(c(0, amplification,
      if(J>2) rep(amplification/2, J-2)), H), H, J, byrow = TRUE,
      dimnames = list(NULL, donors))
  w0 <- W[, T]; wm1 <- W[, T-1]
  beta <- F %*% W
  q_true <- if (scenario %in% c("affine", "hull_exit", "curved")) 2L else 1L
  tail_beta <- t(beta[, (T-q_true+1):(T+H), drop = FALSE])
  energy <- sum(diff(tail_beta, differences = q_true)^2)
  radius <- sqrt(energy)
  z <- trsc_counterfactual_set(population_geometry, deployment, w0, q_true,
       radius, previous = if(q_true==2) wm1 else NULL, type = "population")
  support <- trsc_support(z, rep(1, H))
  # Truth supports are numerical enclosures unless a separate analytic formula
  # is supplied; they are not silently treated as exact endpoint values.
  list(id = scenario, class = if(scenario %in% c("curved", "abrupt"))
         "identification_or_stress_depending_on_declared_R" else "identification_inference",
       parameters = list(seed=seed,T=T,J=J,H=H,rank=r,noise_sd=noise_sd,gap=gap,
         velocity=velocity,amplification=amplification,near_distance=near_distance,
         dependence=dependence),
       observed = list(panel=panel,donors=donors,T=T,x_post=deployment,
         roles=list(calibration=seq_len(floor(T/4)),
                    fitting=(floor(T/4)+3):(T-gap),
                    validation=(T-gap+1):(T-1))),
       truth = list(Gamma_true=crossprod(F),geometry_true=population_geometry,
         beta_true=beta,beta0_true=beta[,T],v0_true=beta[,T]-beta[,T-1],
         w_true=W,anchor_true=w0,previous_true=wm1,x_post=deployment,
         muN_true=rowSums(deployment*t(W[,(T+1):(T+H),drop=FALSE])),
         identified_set_true=z,support_true=support,
         radius_true=radius,q_true=q_true,mu_pre_true=mu_pre))
}

trsc_sim_post_protocol <- function(q, R, predictive_q = c(1L,2L),
    kappa = c(.1,1), bandwidths = 1L, diameter_bound = 2.000001) {
  # Known design bound, fixed before draws. Never replace with an observed
  # oracle error. Broad envelopes intentionally test honest conservative use.
  if (!q %in% c(1,2) || R < 0 || diameter_bound <= 0)
    stop("Invalid maintained protocol.")
  list(q=q,R=R,predictive_q=predictive_q,kappa=kappa,bandwidths=bandwidths,
       diameter_bound=diameter_bound,
       source="Predeclared population profile-diameter bound for bounded-response DGP")
}

trsc_sim_post_feasible <- function(observed, protocol) {
  # Deliberately no scenario/truth argument.
  start <- proc.time()[[3]]
  d <- trsc_data(observed$panel,"unit","period","outcome","Target",
                observed$donors,observed$T+1,seq_len(observed$T),"outcome")
  roles <- observed$roles
  plan <- trsc_sample_plan(d,calibration_periods=roles$calibration,
      fitting_periods=roles$fitting,validation_periods=roles$validation,
      target_period=observed$T)
  workflow <- trsc_workflow(d,plan,q=protocol$predictive_q,kappa=protocol$kappa,
                           bandwidths=protocol$bandwidths)
  g <- workflow$geometry
  if (is.null(g)) g <- workflow$selected$geometry
  # Research adapter: explicit simplex projection for a recorded center.
  # This never changes Estimator A. Uniform diameter envelopes below apply to
  # the projected center itself; no fitted-risk claim is made for this adapter.
  normalize <- function(w) {
    w <- pmax(w,0); if(sum(w)==0) stop("Invalid fitted representative.")
    w/sum(w)
  }
  candidates <- workflow$candidate_library
  fits <- attr(candidates,"fits")
  ids <- names(candidates)
  D <- protocol$diameter_bound
  rho <- if(protocol$q==1) D else sqrt(5)*D
  certificates <- lapply(ids,function(id) {
    readout <- candidates[[id]]
    fit <- fits[[readout$candidate_specification$fit_name]]
    center <- normalize(readout$representative)
    previous <- normalize(fit$weights[,ncol(fit$weights)-1])
    trsc_anchor_certificate(g,center,protocol$q,
       c(estimation=rho,fitting_gap=rho,causal_gap=0,numerical=0),0,
       protocol$source,previous=if(protocol$q==2) previous else NULL,
       provenance="maintained_bound",event_id="bounded_design_all_candidates",
       candidate_id=id,covered_candidates=ids)
  })
  names(certificates) <- ids
  ac <- certificates[[workflow$selection$selected_name]]
  # Absolute all-chord norm difference <= sum of the two chord diameters.
  # Empirical diameter is bounded from observable column norms. A positive
  # numerical cushion prevents a research-level rounded inequality from being
  # advertised as an exact arithmetic bound.
  empirical_bound <- 2*max(sqrt(colSums(g$B^2)))
  epsilon <- (D+empirical_bound)*(1+1e-10)
  gc <- trsc_geometry_certificate(g,epsilon,0,protocol$source,
                                 provenance="maintained_bound")
  center <- as.numeric(ac$anchors[[length(ac$anchors)]]);names(center)<-observed$donors
  previous <- if(protocol$q==2) as.numeric(ac$anchors[[1]]) else NULL
  outer <- trsc_confidence_set(observed$x_post,ac,gc,protocol$q,protocol$R,
                              library=certificates)
  plugin <- trsc_counterfactual_set(g,observed$x_post,center,protocol$q,protocol$R,
                                   previous=previous)
  list(layer="feasible_maintained_bounds",workflow=workflow,center=center,
       previous=previous,plugin=plugin,outer=outer,protocol=protocol,
       candidate_certificates=certificates,geometry_certificate=gc,
       runtime=proc.time()[[3]]-start,
       center_convention="explicit simplex projection; supplied envelopes apply to recorded centers")
}

trsc_sim_post_oracle <- function(scenario, R, q = scenario$truth$q_true) {
  trsc_counterfactual_set(scenario$truth$geometry_true,scenario$observed$x_post,
      scenario$truth$anchor_true,q,R,
      previous=if(q==2) scenario$truth$previous_true else NULL,type="population")
}

trsc_sim_post_semifeasible <- function(run, truth) {
  # Explicit oracle geometry; fitted centers and predeclared diameter envelopes.
  old <- run$outer$anchor_certificate
  ac <- trsc_anchor_certificate(truth$geometry_true,run$center,run$protocol$q,
      old$components,0,"Semi-feasible: known population geometry and maintained diameter",
      previous=run$previous)
  gc <- trsc_geometry_certificate(truth$geometry_true,0,0,
      "Oracle population geometry",provenance="analytically_known")
  out <- trsc_confidence_set(truth$x_post,ac,gc,run$protocol$q,run$protocol$R)
  attr(out,"simulation_layer") <- "semi_feasible_oracle_geometry"
  out
}

trsc_sim_post_metrics <- function(run, truth, a = rep(1,nrow(truth$x_post))) {
  started <- proc.time()[[3]]
  pop <- trsc_counterfactual_set(truth$geometry_true,truth$x_post,
    truth$anchor_true,run$protocol$q,run$protocol$R,
    previous=if(run$protocol$q==2) truth$previous_true else NULL,type="population")
  hi <- trsc_support(run$outer,a);lo <- trsc_support(run$outer,-a)
  phi <- trsc_support(pop,a);plo <- trsc_support(pop,-a)
  bands <- trsc_bounds(run$outer)
  coordinate_pop <- trsc_bounds(pop)
  ac <- run$outer$anchor_certificate;F <- truth$geometry_true$B
  e0 <- drop(F %*% (run$center-truth$anchor_true))
  ev <- if(run$protocol$q==2) drop(F %*%
    (run$center-run$previous-truth$anchor_true+truth$previous_true)) else 0
  joint <- sqrt(sum(e0^2)+sum(ev^2))
  anchor_good <- if(run$protocol$q==1) sqrt(sum(e0^2))<=ac$radius0 else joint<=ac$joint_radius
  # Conservative sufficient geometry event, evaluated using truth only here.
  chord_sum <- 2*max(sqrt(colSums(F^2))) + 2*max(sqrt(colSums(run$outer$geometry$B^2)))
  geometry_good <- chord_sum <= run$outer$geometry_certificate$epsilon
  whole_event <- anchor_good && geometry_good
  # Exact support comparison is tri-state if oracle optimizer is not exact.
  comparison <- function(upper,oracle) {
    if(upper >= oracle$upper) return("certified_covered")
    if(is.finite(oracle$primal_lower) && upper < oracle$primal_lower)
      return("certified_not_covered")
    "unresolved"
  }
  rmin <- trsc_min_radius(truth$geometry_true,truth$anchor_true,
    run$protocol$q,nrow(truth$x_post),
    previous=if(run$protocol$q==2) truth$previous_true else NULL)
  hull <- sum(abs(a)*apply(truth$x_post,1,function(x)diff(range(x))))
  infeasible <- phi$status == "structural_infeasible" || plo$status == "structural_infeasible"
  pop_width <- if(infeasible) NA_real_ else phi$upper+plo$upper
  confidence_width <- hi$upper+lo$upper
  q <- run$protocol$q;H <- nrow(truth$x_post)
  beta_tail <- t(tail(t(truth$beta_true),H+q))
  actual_energy <- sum(diff(t(beta_tail),differences=q)^2)
  data.frame(layer=run$layer,cutoff_profile_error=sum(e0^2),
    structural_infeasible=infeasible,
    continuation_in_class=actual_energy<=run$protocol$R^2+1e-14,
    population_width_enclosure=pop_width,confidence_width=confidence_width,
    numerical_enclosure_width_difference=confidence_width-pop_width,
    upper_support_coverage=comparison(hi$upper,phi),
    lower_support_coverage=comparison(lo$upper,plo),
    true_path_in_coordinate_band=all(truth$muN_true>=bands$lower & truth$muN_true<=bands$upper),
    whole_set_sufficient_event=whole_event,whole_set_coverage_exact=NA,
    Rmin_lower=rmin$lower,Rmin_upper=rmin$upper,
    datewise_relaxed_upper=sum(ifelse(a>=0,a*coordinate_pop$upper,a*coordinate_pop$lower)),
    joint_population_upper_enclosure=phi$upper,
    width_to_hull=if(hull>0)confidence_width/hull else NA_real_,
    certified_upper=hi$certificate$valid && lo$certificate$valid,
    certified_gap=hi$certified_gap,solver_status=hi$status,
    certificate_method=hi$certificate$method,
    runtime=run$runtime+proc.time()[[3]]-started,stringsAsFactors=FALSE)
}

trsc_sim_post_grid <- function(scenarios, seeds, protocol, ...) {
  metrics <- list(); failures <- list(); attempts <- list(); k <- 0L
  for(id in scenarios) for(seed in seeds) {
    k <- k+1L
    warnings <- character()
    result <- tryCatch(withCallingHandlers({
      d <- trsc_sim_post_scenario(id,seed=seed,...)
      run <- trsc_sim_post_feasible(d$observed,protocol)
      trsc_sim_post_metrics(run,d$truth)
    }, warning=function(w) { warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning") }),
      error=function(e)e)
    attempts[[k]] <- data.frame(scenario=id,seed=seed,warnings=paste(warnings,collapse=" | "))
    if(inherits(result,"error")) failures[[length(failures)+1L]] <-
      data.frame(scenario=id,seed=seed,error_class=class(result)[1],message=conditionMessage(result))
    else metrics[[length(metrics)+1L]] <- cbind(scenario=id,seed=seed,result)
  }
  bind <- function(x) if(length(x))do.call(rbind,x) else data.frame()
  git_sha <- tryCatch(system2("git",c("rev-parse","HEAD"),stdout=TRUE,stderr=FALSE),error=function(e)NA_character_)
  list(metrics=bind(metrics),failures=bind(failures),attempts=bind(attempts),
       manifest=list(attempted=k,successful=length(metrics),failed=length(failures),
         protocol=protocol,scenarios=scenarios,seeds=seeds,parameters=list(...),
         git_sha=git_sha,dirty=length(system2("git",c("status","--porcelain"),stdout=TRUE))>0,
         R=R.version.string,package=as.character(utils::packageVersion("TRSC")),
         versions=vapply(c("osqp","ECOSolveR","gmp"),function(p)as.character(utils::packageVersion(p)),character(1)),
         RNGkind=RNGkind(),session=utils::sessionInfo(),scientific_findings_claimed=FALSE))
}

trsc_sim_post_report <- function(grid) {
  # Preserve all attempts. Binomial summaries concern MC repetition frequency,
  # not an extra confidence statement about a causal estimand.
  stopifnot(grid$manifest$attempted == nrow(grid$metrics)+nrow(grid$failures))
  x <- grid$metrics
  numeric_metrics <- c("cutoff_profile_error","population_width_enclosure",
                       "confidence_width","numerical_enclosure_width_difference",
                       "width_to_hull","certified_gap","runtime")
  summaries <- list();coverage <- list()
  if(nrow(x)) for(id in unique(x$scenario)) {
    d <- x[x$scenario==id,,drop=FALSE]
    for(m in numeric_metrics) {
      v <- d[[m]];v <- v[is.finite(v)]
      summaries[[length(summaries)+1L]] <- data.frame(scenario=id,metric=m,n=length(v),
        mean=if(length(v))mean(v) else NA_real_,
        mcse=if(length(v)>1)sd(v)/sqrt(length(v)) else NA_real_)
    }
    for(m in c("true_path_in_coordinate_band","whole_set_sufficient_event")) {
      v <- d[[m]];n <- sum(!is.na(v));k <- sum(v,na.rm=TRUE)
      ci <- if(n) stats::binom.test(k,n)$conf.int else c(NA_real_,NA_real_)
      coverage[[length(coverage)+1L]] <- data.frame(scenario=id,metric=m,
        successes=k,n=n,frequency=if(n)k/n else NA_real_,
        mc_lower=ci[1],mc_upper=ci[2],
        interpretation=if(m=="whole_set_sufficient_event")
          "lower bound from sufficient events, not exact whole-set coverage" else
          "simultaneous coordinate-band coverage of the true mean path")
    }
  }
  bind <- function(z) if(length(z))do.call(rbind,z) else data.frame()
  list(summary=bind(summaries),coverage=bind(coverage),failures=grid$failures,
       attempts=grid$attempts,manifest=grid$manifest,
       unresolved_endpoints=if(nrow(x))sum(x$upper_support_coverage=="unresolved"|
                                          x$lower_support_coverage=="unresolved") else 0L)
}

trsc_sim_post_write <- function(report,directory) {
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  for(n in c("summary","coverage","failures","attempts"))
    utils::write.csv(report[[n]],file.path(directory,paste0(n,".csv")),row.names=FALSE)
  saveRDS(report$manifest,file.path(directory,"manifest.rds"))
  writeLines(capture.output(report$manifest$session),file.path(directory,"sessionInfo.txt"))
  invisible(directory)
}

library(TRSC)
source("research/simulations/R/scenarios.R")
source("research/simulations/R/metrics.R")
source("research/simulations/R/benchmarks.R")
source("research/simulations/R/run-study.R")
source("research/simulations/R/post-treatment.R")
source("research/simulations/R/frozen-study.R")
f<-trsc_sim_validate_frozen_registries()
stopifnot(f$cell_count==179L,f$attempted_repetitions==73400L,f$registry_hash=="89181ec928012f3424e5980d540b8af1d89a23be83afe9d2edf1ae9e367b5b91")
h1<-trsc_sim_registry_hashes();h2<-trsc_sim_registry_hashes();stopifnot(identical(h1,h2))
sizes<-table(f$candidates$library_id);stopifnot(identical(as.integer(sizes[c("L1","L4","L12","L24")]),c(1L,4L,12L,24L)))
stopifnot(identical(as.numeric(f$epsilon$epsilon_over_D),c(0,.01,.1,.5,1)))
l12<-trsc_sim_candidate_spec("L12",f);stopifnot(identical(l12$q,c(1L,2L)),identical(l12$kappa,c(.03,.3,3)),identical(l12$bandwidths,c(1L,4L)))
cells<-f$cells;stopifnot(!anyDuplicated(cells$cell_id),identical(as.integer(cells$cell_index),seq_len(nrow(cells))))
stopifnot(trsc_sim_frozen_seed(14,1)==1140001L)
a<-trsc_sim_attempt_table(f);stopifnot(nrow(a)==73400L,!anyDuplicated(a$primary_key),!any(a$seed%in%c(11,12,101)))
stopifnot(any(cells$panel=="E1"&cells$J_requested==120),any(cells$panel=="P2"&cells$near_distance==0),any(cells$panel=="P4"&cells$radius_rule=="zero"))
fr<-cells[cells$panel=="P0"&cells$scenario=="full_rank"&cells$structural_q==1&cells$radius_rule=="rmin",]
stopifnot(nrow(fr)==1,fr$rank_requested==3,fr$rank_effective==10,grepl("forces effective rank",fr$generator_override))
d<-trsc_sim_post_scenario("full_rank",seed=11,T=60,J=10,H=2,rank=3);stopifnot(identical(.trsc_verify_dimensions(fr,d),list(J=10L,rank=10L)))
fc<-cells[cells$panel=="E0"&cells$scenario=="face_change",];stopifnot(fc$J_requested==10,fc$J_effective==4)
dx<-trsc_sim_scenario("face_change",T=60,J=10,seed=11,reference_size=100);stopifnot(dx$J==4L)
ed<-cells[cells$panel=="E0"&cells$scenario=="exact_duplication",];stopifnot(ed$J_effective==11)
stopifnot(identical(names(formals(trsc_sim_post_feasible)),c("observed","protocol")))
dummy<-list(geometry_true=trsc_geometry(matrix(c(0,1),1,dimnames=list(NULL,c("D1","D2")))),anchor_true=c(D1=.5,D2=.5),previous_true=c(D1=.5,D2=.5))
r<-cells[cells$panel=="P0",][1,];r$structural_q<-1;r$H<-5;r$D<-2.000001;r$radius_rule<-"rmin_plus_D"
rr<-trsc_sim_resolve_radius(r,dummy);stopifnot(rr$Rmin_lower==0,rr$Rmin_upper==0,abs(rr$R-2.000001)<1e-12)
dh<-trsc_sim_post_scenario("hull_exit",seed=12,T=60,J=3,H=2,rank=1)
r2<-cells[cells$panel=="P0"&cells$scenario=="hull_exit"&cells$structural_q==2,][1,];r2$H<-2;r2$radius_rule<-"zero"
rr2<-trsc_sim_resolve_radius(r2,dh$truth);stopifnot(rr2$R==0,rr2$Rmin_upper>0);if(rr2$Rmin_lower>0)stopifnot(rr2$structural_infeasible)
dirs<-trsc_sim_support_directions(3);stopifnot(all(c("e1","minus_e1","e2","minus_e2","e3","minus_e3","all_ones","minus_all_ones","equal_weight_mean","minus_equal_weight_mean")%in%names(dirs)))
pb<-trsc_sim_probability_budget(.05,4,.01,.005,.005,.001);stopifnot(abs(pb$alpha_anchor-.03)<1e-12,abs(pb$per_candidate_anchor_budget-.0075)<1e-12,abs(pb$alpha_path-.0065)<1e-12)
bad<-tryCatch({trsc_sim_validate_nominal_risk_inputs(list(path_risk_bound=1));NULL},error=identity);stopifnot(inherits(bad,"error"),grepl("restricted_coercivity_bound",conditionMessage(bad)))
stopifnot(trsc_sim_nominal_inference_gate()$status=="NOT READY",trsc_sim_execution_gate(f)$status=="READY")
tmp<-tempfile("combine-");dir.create(tmp);small<-f;small$cells<-cells[1,,drop=FALSE];small$cells$repetitions<-2L;ex<-trsc_sim_attempt_table(small)
a1<-trsc_sim_frozen_attempt(small$cells,1,small);a2<-trsc_sim_frozen_attempt(small$cells,2,small)
for(x in list(a1,a2)){x$status<-if(x$repetition==1)"success"else"error";x$error_class<-if(x$repetition==1)NA_character_ else "fixture";x$message<-NA_character_;x$warnings<-"";x$runtime_seconds<-0;x$rng_kind<-"fixture";trsc_sim_atomic_save(list(attempt=x,result=NULL),trsc_sim_result_path(tmp,x))}
z<-trsc_sim_combine_frozen(tmp,small,ex);stopifnot(nrow(z$attempts)==2,sum(z$attempts$status=="error")==1)
unlink(trsc_sim_result_path(tmp,a2));e<-tryCatch({trsc_sim_combine_frozen(tmp,small,ex);NULL},error=identity);stopifnot(inherits(e,"error"),grepl("missing",conditionMessage(e),ignore.case=TRUE));unlink(tmp,recursive=TRUE)


# Pre-run provenance must be distinct from post-build workspace dirtiness.
prov<-trsc_sim_provenance_snapshot(initial_sha="0123456789abcdef",initial_dirty=FALSE)
stopifnot(prov$repository_sha=="0123456789abcdef",identical(prov$dirty,FALSE),
          all(c("post_check_dirty","post_check_dirty_paths")%in%names(prov)))

# Reserved-seed reporting fixture: require the complete predeclared numerical
# diagnostics without executing a scientific registry repetition.
pc<-cells[cells$panel=="P0"&cells$scenario=="persistence"&cells$structural_q==1&
          cells$radius_rule=="rmin",][1,,drop=FALSE]
pc$T<-60;pc$J_requested<-3;pc$J_effective<-3;pc$rank_requested<-1;pc$rank_effective<-1
pc$H<-1;pc$gap<-10;pc$velocity_rule<-"0.05/T";pc$velocity<-.05/60
pc$predictive_library<-"L1"
pd<-trsc_sim_post_scenario("persistence",seed=12,T=60,J=3,H=1,rank=1,
                           noise_sd=.2,gap=10,velocity=.05/60,
                           amplification=1,near_distance=.001,dependence=0)
prad<-trsc_sim_resolve_radius(pc,pd$truth)
pp<-trsc_sim_post_protocol(1,prad$R,predictive_q=1,kappa=.3,bandwidths=1,
                           diameter_bound=pc$D)
pfit<-trsc_sim_post_feasible(pd$observed,pp)
pm<-trsc_sim_post_frozen_metrics(pc,pfit,pd$truth,prad,f)
required_support<-c("upper_dual_objective","lower_dual_objective",
 "upper_certified_correction","lower_certified_correction",
 "upper_dual_feasible","lower_dual_feasible",
 "unrestricted_hull_lower","unrestricted_hull_upper",
 "unrestricted_hull_width","width_to_unrestricted_hull",
 "lower_endpoint_comparison","upper_endpoint_comparison")
stopifnot(all(required_support%in%names(pm$supports)))
allowed<-c("certified_covered","certified_not_covered","unresolved")
stopifnot(all(pm$supports$lower_endpoint_comparison%in%allowed),
          all(pm$supports$upper_endpoint_comparison%in%allowed))
required_decomp<-c("direction","population_width","plugin_width",
                   "semi_feasible_width","outer_width","statistical_enlargement")
stopifnot(is.data.frame(pm$decomposition),
          all(required_decomp%in%names(pm$decomposition)))

cat("Frozen Stage-A registry and orchestration contracts passed.\n")

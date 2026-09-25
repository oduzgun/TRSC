library(TRSC)
source("research/simulations/R/post-treatment.R")

for (id in c("persistence", "affine", "hull_exit", "exact_ambiguity",
             "near_redundancy", "weak_excitation", "full_rank", "curved", "abrupt")) {
  d <- trsc_sim_post_scenario(id, seed = 101, T = 60, J = 3, H = 2)
  stopifnot(all(c("Gamma_true", "beta_true", "beta0_true", "v0_true", "w_true",
                 "x_post", "muN_true", "identified_set_true", "support_true") %in% names(d$truth)))
  oracle <- trsc_sim_post_oracle(d, R = d$truth$radius_true)
  stopifnot(inherits(oracle, "trsc_identified_set"))
}

d <- trsc_sim_post_scenario("persistence", seed = 11, T = 60, J = 3, H = 2)
p <- trsc_sim_post_protocol(q = 1, R = 0)
run <- trsc_sim_post_feasible(d$observed, p)
stopifnot(run$layer == "feasible_maintained_bounds",
          run$outer$coverage$failure_probability == 0,
          isTRUE(trsc_sim_post_metrics(run, d$truth)$whole_set_sufficient_event))
# Mutation of truth cannot affect the observed-input runner.
altered <- d; altered$truth$muN_true[] <- 1e9; altered$truth$Gamma_true[] <- -999
repeat_run <- trsc_sim_post_feasible(altered$observed, p)
stopifnot(identical(run$center, repeat_run$center),
          identical(run$outer$energy_allowance, repeat_run$outer$energy_allowance),
          identical(trsc_bounds(run$outer)[1:3], trsc_bounds(repeat_run$outer)[1:3]))

for (id in c("affine", "hull_exit", "exact_ambiguity")) {
  d <- trsc_sim_post_scenario(id, seed = 12, T = 60, J = 2, H = 2,
                             amplification = 100)
  q <- if (id == "exact_ambiguity") 1 else 2
  p <- trsc_sim_post_protocol(q = q, R = if (id == "hull_exit") 1 else 0)
  r <- trsc_sim_post_feasible(d$observed, p)
  metrics <- trsc_sim_post_metrics(r, d$truth)
  stopifnot(metrics$true_path_in_coordinate_band,
            metrics$whole_set_sufficient_event,
            is.finite(metrics$confidence_width),
            metrics$certified_upper,
            is.na(metrics$whole_set_coverage_exact))
}

grid <- trsc_sim_post_grid(scenarios = c("persistence", "not_a_design"),
                          seeds = 11, T = 60, J = 2, H = 1,
                          protocol = trsc_sim_post_protocol(1, 0))
stopifnot(grid$manifest$attempted == 2, nrow(grid$failures) == 1,
          nrow(grid$metrics) == 1, !grid$manifest$scientific_findings_claimed)
report <- trsc_sim_post_report(grid)
stopifnot(nrow(report$coverage)==2, all(report$coverage$n==1))
folder <- tempfile("post-smoke-")
trsc_sim_post_write(report,folder)
stopifnot(file.exists(file.path(folder,"manifest.rds")),
          file.exists(file.path(folder,"failures.csv")))
unlink(folder,recursive=TRUE)
cat("Post-treatment simulation smoke checks passed.\n")

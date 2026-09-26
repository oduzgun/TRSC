library(TRSC)
source("research/simulations/R/scenarios.R")
source("research/simulations/R/metrics.R")
source("research/simulations/R/benchmarks.R")
source("research/simulations/R/run-study.R")
source("research/simulations/R/post-treatment.R")
source("research/simulations/R/frozen-study.R")

frozen <- trsc_sim_validate_frozen_registries()
stopifnot(
  frozen$cell_count == 179L,
  frozen$attempted_repetitions == 73400,
  identical(frozen$registry_hash,
            "fb0d15a82c254695cc15e7709780981bfc405d3b8ccb6e34d3c2f347d821eece")
)

sizes <- table(frozen$candidates$library_id)
stopifnot(identical(as.integer(sizes[c("L1", "L4", "L12", "L24")]),
                    c(1L, 4L, 12L, 24L)))
stopifnot(identical(as.numeric(frozen$epsilon$epsilon_over_D), c(0, .01, .1, .5, 1)))

cells <- frozen$cells
stopifnot(!anyDuplicated(cells$cell_id),
          identical(cells$cell_index, seq_len(nrow(cells))))

stopifnot(trsc_sim_frozen_seed(14L, 1L) == 1140001L)
key <- trsc_sim_attempt_key(frozen$registry_hash, cells$cell_id[1], 1L)
stopifnot(identical(key, paste(frozen$registry_hash, cells$cell_id[1], 1L, sep = ":")))

full_rank <- cells[cells$panel == "P0" & cells$scenario == "full_rank" &
                   cells$structural_q == 1 & cells$radius_rule == "rmin", ]
stopifnot(nrow(full_rank) == 1L,
          full_rank$rank_requested == 3,
          full_rank$rank_effective == 10,
          grepl("forces effective rank", full_rank$generator_override))
face_change <- cells[cells$panel == "E0" & cells$scenario == "face_change", ]
stopifnot(nrow(face_change) == 1L, face_change$J_requested == 10,
          face_change$J_effective == 4)

p0 <- cells[cells$panel == "P0", ]
stopifnot(all(p0$structural_q %in% c(1, 2)), all(p0$predictive_library == "L12"))
l12 <- trsc_sim_candidate_spec("L12", frozen)
stopifnot(identical(l12$q, c(1L, 2L)),
          identical(l12$kappa, c(.03, .3, 3)),
          identical(l12$bandwidths, c(1L, 4L)))

dummy_truth <- list(
  geometry_true = trsc_geometry(matrix(c(0, 1), 1,
    dimnames = list(NULL, c("D1", "D2")))),
  anchor_true = c(D1 = .5, D2 = .5),
  previous_true = c(D1 = .5, D2 = .5)
)
r <- p0[1, ]; r$structural_q <- 1; r$H <- 5; r$D <- 2.000001; r$radius_rule <- "rmin_plus_D"
rr <- trsc_sim_resolve_radius(r, dummy_truth)
stopifnot(rr$Rmin_lower == 0, rr$Rmin_upper == 0,
          abs(rr$R - 2.000001) < 1e-12,
          grepl("not feasible data-driven", rr$provenance))

dirs <- trsc_sim_support_directions(3)
stopifnot(all(c("e1", "minus_e1", "e2", "minus_e2", "e3", "minus_e3",
                "all_ones", "minus_all_ones", "equal_weight_mean",
                "minus_equal_weight_mean") %in% names(dirs)))

pb <- trsc_sim_probability_budget(alpha = .05, M = 4, alpha_geometry = .01,
                                  alpha_mixing = .005, alpha_other = .005,
                                  fitting_failure = .001)
stopifnot(abs(pb$alpha_anchor - .03) < 1e-12,
          abs(pb$per_candidate_anchor_budget - .0075) < 1e-12,
          abs(pb$alpha_path - .0065) < 1e-12)

ig <- trsc_sim_nominal_inference_gate()
stopifnot(identical(ig$status, "NOT READY"),
          "restricted_coercivity_bound" %in% ig$missing_inputs,
          "fitting_to_cutoff_bound" %in% ig$missing_inputs)
eg <- trsc_sim_execution_gate(frozen)
stopifnot(identical(eg$status, "READY"))

attempts <- trsc_sim_attempt_table(frozen)
stopifnot(nrow(attempts) == 73400L,
          !anyDuplicated(attempts$primary_key),
          !any(attempts$seed %in% c(11L, 12L, 101L)))

tmp <- tempfile("trsc-frozen-combine-")
dir.create(tmp)
small <- frozen
small$cells <- frozen$cells[1, , drop = FALSE]
small$cells$repetitions <- 2L
small$registry_hash <- frozen$registry_hash
fake1 <- trsc_sim_frozen_attempt(small$cells[1, ], 1L, small)
trsc_sim_atomic_save(list(attempt = c(fake1, list(status = "error")), result = NULL),
                     trsc_sim_result_path(tmp, fake1))
err <- tryCatch({ trsc_sim_combine_frozen(tmp, small); NULL }, error = identity)
stopifnot(inherits(err, "error"), grepl("missing", conditionMessage(err), ignore.case = TRUE))
unlink(tmp, recursive = TRUE)

cat("Frozen Stage-A registry and orchestration contracts passed.\n")

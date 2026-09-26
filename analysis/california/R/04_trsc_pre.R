source(file.path("analysis", "california", "R", "00_helpers.R"))
library(TRSC)

panel <- read_processed_long()
panel_pre <- panel[panel$year <= 1988, ]
donors <- sort(setdiff(unique(panel_pre$state), "California"))

prepared <- trsc_data(
  panel_pre, unit = "state", time = "year", outcome = "cigsale",
  treated = "California", donors = donors,
  intervention_time = 1989, pre_period = 1970:1988,
  predictors = "cigsale", period_step = 1
)
plan <- trsc_sample_plan(
  prepared,
  calibration_periods = 1970:1974,
  fitting_periods = 1975:1985,
  validation_periods = 1986:1988,
  calibration_thin = 1, validation_thin = 1,
  target_period = 1988, min_gap0 = 1, min_gap1 = 1
)
geometry <- trsc_calibration_geometry(prepared, plan)
saveRDS(prepared, file.path(app_root, "data", "processed", "trsc_prepared.rds"))
saveRDS(plan, file.path(app_root, "data", "processed", "trsc_plan.rds"))
saveRDS(geometry, file.path(app_root, "data", "processed", "calibration_geometry.rds"))

initial_kappa <- c(.03, .1, .3, 1, 3, 10, 30)
bandwidths <- c(1L, 2L, 4L)

run_order <- function(q_value) {
  grid <- initial_kappa
  for (attempt in seq_len(10L)) {
    candidates <- trsc_candidate_library(
      prepared, plan, kappa = grid, q = q_value,
      bandwidths = bandwidths, geometry = geometry
    )
    selection <- trsc_select(
      prepared, candidates,
      validation_periods = 1986:1988,
      target_period = 1988
    )
    selected_kappa <- selection$selected$candidate_specification$kappa
    if (selected_kappa > min(grid) && selected_kappa < max(grid)) break
    if (selected_kappa == min(grid)) grid <- sort(unique(c(min(grid) * .1, grid)))
    if (selected_kappa == max(grid)) grid <- sort(unique(c(grid, max(grid) * 10)))
    if (attempt == 10L) stop("Kappa boundary expansion did not select an interior value.")
  }
  fit <- selected_fit(selection)
  spec <- selection$selected$candidate_specification
  scores <- data.frame(
    candidate = names(selection$scores), score = as.numeric(selection$scores),
    q = q_value, stringsAsFactors = FALSE
  )
  grid_info <- attr(selection$candidates, "grid")
  scores <- merge(scores, grid_info, by = c("candidate", "q"), all.x = TRUE, sort = FALSE)
  write.csv(scores, file.path(output_root, "diagnostics", paste0("q", q_value, "_candidate_scores.csv")), row.names = FALSE)
  write.csv(data.frame(
    q = q_value, selected_candidate = selection$selected_name,
    kappa = spec$kappa, h = spec$h,
    validation_score = selection$scores[[selection$selected_name]],
    fit_rmse = sqrt(mean(fit$residuals^2)),
    fit_loss = fit$objective$fit_loss,
    temporal_roughness = fit$objective$temporal_roughness,
    penalty = fit$objective$regularization_penalty,
    fitting_cutoff = max(fit$periods),
    solver_status = fit$solver$status,
    iterations = fit$solver$iterations,
    primal_residual = fit$solver$primal_residual,
    dual_residual = fit$solver$dual_residual,
    stringsAsFactors = FALSE
  ), file.path(output_root, "tables", paste0("q", q_value, "_selection.csv")), row.names = FALSE)
  saveRDS(selection, file.path(output_root, "diagnostics", paste0("q", q_value, "_selection.rds")))
  saveRDS(fit, file.path(output_root, "diagnostics", paste0("q", q_value, "_selected_fit.rds")))
  write.csv(data.frame(period = fit$periods, fitted = fit$fitted_values,
                       observed = fit$target_outcomes, residual = fit$residuals),
            file.path(output_root, "diagnostics", paste0("q", q_value, "_fitted_path.csv")), row.names = FALSE)
  write.csv(as.data.frame(fit$weights),
            file.path(output_root, "diagnostics", paste0("q", q_value, "_representative_weights.csv")))
  write.csv(as.data.frame(fit$metric_coordinates),
            file.path(output_root, "diagnostics", paste0("q", q_value, "_metric_coordinates.csv")))
  invisible(selection)
}

q1 <- run_order(1L)
q2 <- run_order(2L)
write.csv(data.frame(year = 1986:1988, validation_weight = q1$validation_weights),
          file.path(output_root, "diagnostics", "validation_weights.csv"), row.names = FALSE)
write_status("04_trsc_pre", "completed")
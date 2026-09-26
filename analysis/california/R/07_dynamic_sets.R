source(file.path("analysis", "california", "R", "00_helpers.R"))
library(TRSC)

wide <- read.csv(processed_wide_path, check.names = FALSE)
donors <- setdiff(names(wide), c("year", "California"))
post <- wide[wide$year >= 1989, ]
Xpost <- as.matrix(post[donors])
colnames(Xpost) <- donors
H <- nrow(Xpost)
geometry <- readRDS(file.path(app_root, "data", "processed", "calibration_geometry.rds"))
q1_anchor <- readRDS(file.path(output_root, "diagnostics", "q1_boundary_level.rds"))
q2_boundary <- readRDS(file.path(output_root, "diagnostics", "q2_boundary.rds"))
q1_fit <- readRDS(file.path(output_root, "diagnostics", "q1_selected_fit.rds"))
q2_fit <- readRDS(file.path(output_root, "diagnostics", "q2_selected_fit.rds"))

increment_scale <- function(fit, q_value, H) {
  coordinates <- t(fit$metric_coordinates)
  increments <- coordinates[-1, , drop = FALSE] - coordinates[-nrow(coordinates), , drop = FALSE]
  if (q_value == 2L) {
    increments <- increments[-1, , drop = FALSE] - increments[-nrow(increments), , drop = FALSE]
  }
  root_energy <- sqrt(rowSums(increments^2))
  list(scale = sqrt(H) * median(root_energy), increments = root_energy)
}

hull <- data.frame(
  year = post$year,
  lower = apply(Xpost, 1, min),
  upper = apply(Xpost, 1, max)
)
hull$width <- hull$upper - hull$lower
write.csv(hull, file.path(output_root, "diagnostics", "unrestricted_donor_hull.csv"), row.names = FALSE)

run_order <- function(q_value, anchor, previous, fit) {
  rmin <- trsc_min_radius(
    geometry, anchor = anchor, q = q_value, H = H, previous = previous
  )
  scale_info <- increment_scale(fit, q_value, H)
  c_values <- c(0, .25, .5, 1, 2, 4, 8, 16)
  if (scale_info$scale == 0) c_values <- 0

  rows <- list(); logs <- list(); path_rows <- list(); counter <- 0L
  saturated <- FALSE
  repeat {
    for (c_value in setdiff(c_values, vapply(rows, function(x) x$c[[1]], numeric(1)))) {
      R <- rmin$upper + c_value * scale_info$scale
      object <- trsc_counterfactual_set(
        geometry, Xpost, anchor = anchor, q = q_value, R = R,
        previous = previous, type = "plugin"
      )
      bounds <- trsc_bounds(object)
      bound_width <- bounds$upper - bounds$lower
      results <- attr(bounds, "support_results")
      statuses <- unlist(lapply(results, function(z) c(z$lo$status, z$hi$status)))
      feasible <- !any(statuses == "structural_infeasible")
      rows[[length(rows) + 1L]] <- data.frame(
        q = q_value, c = c_value, R = R,
        horizon = bounds$horizon, year = post$year,
        lower = bounds$lower, upper = bounds$upper, width = bound_width,
        hull_width = hull$width, hull_fraction = bound_width / hull$width,
        feasible = feasible, stringsAsFactors = FALSE
      )
      for (h in seq_len(H)) {
        counter <- counter + 1L
        logs[[counter]] <- support_record(results[[h]]$lo, q_value, R, h, "lower", post$year[h])
        counter <- counter + 1L
        logs[[counter]] <- support_record(results[[h]]$hi, q_value, R, h, "upper", post$year[h])
      }
      average_support <- trsc_support(object, rep(1 / H, H))
      total_support <- trsc_support(object, rep(1, H))
      path_rows[[length(path_rows) + 1L]] <- data.frame(
        q = q_value, c = c_value, R = R,
        average_upper = average_support$upper,
        total_upper = total_support$upper,
        average_status = average_support$status,
        total_status = total_support$status,
        stringsAsFactors = FALSE
      )
    }
    combined <- do.call(rbind, rows)
    last_c <- max(combined$c)
    last <- combined[combined$c == last_c, ]
    saturated <- nrow(last) == H && all(abs(last$width - last$hull_width) <= 1e-6 * pmax(1, last$hull_width))
    if (saturated || scale_info$scale == 0 || max(c_values) >= 1024) break
    c_values <- c(c_values, max(c_values) * 2)
  }

  write.csv(data.frame(q = q_value, lower = rmin$lower, upper = rmin$upper,
                       estimate = rmin$estimate, status = rmin$status,
                       scale = scale_info$scale, saturated = saturated),
            file.path(output_root, "tables", paste0("q", q_value, "_rmin_scale.csv")), row.names = FALSE)
  write.csv(data.frame(q = q_value, increment = seq_along(scale_info$increments),
                       norm = scale_info$increments),
            file.path(output_root, "diagnostics", paste0("q", q_value, "_predictive_increment_norms.csv")), row.names = FALSE)
  write.csv(do.call(rbind, rows),
            file.path(output_root, "tables", paste0("q", q_value, "_dynamic_bounds.csv")), row.names = FALSE)
  write.csv(do.call(rbind, path_rows),
            file.path(output_root, "tables", paste0("q", q_value, "_path_functionals.csv")), row.names = FALSE)
  write.csv(do.call(rbind, logs),
            file.path(output_root, "solver-logs", paste0("q", q_value, "_dynamic_support_solver.csv")), row.names = FALSE)
}

run_order(1L, q1_anchor, NULL, q1_fit)
run_order(2L, q2_boundary$level, q2_boundary$previous, q2_fit)
write_status("07_dynamic_sets", "completed")
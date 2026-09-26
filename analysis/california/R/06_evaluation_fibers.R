source(file.path("analysis", "california", "R", "00_helpers.R"))
library(TRSC)

wide <- read.csv(processed_wide_path, check.names = FALSE)
donors <- setdiff(names(wide), c("year", "California"))
post <- wide[wide$year >= 1989, ]
Xpost <- as.matrix(post[donors])
colnames(Xpost) <- donors

primary_geometry <- readRDS(file.path(app_root, "data", "processed", "calibration_geometry.rds"))
full_pre_geometry <- trsc_geometry(as.matrix(wide[wide$year <= 1988, donors, drop = FALSE]))
saveRDS(full_pre_geometry, file.path(output_root, "diagnostics", "full_pre_geometry_diagnostic.rds"))

q1_anchor <- readRDS(file.path(output_root, "diagnostics", "q1_boundary_level.rds"))
q2_boundary <- readRDS(file.path(output_root, "diagnostics", "q2_boundary.rds"))

run_fibers <- function(geometry, geometry_role, anchor, q_value) {
  rows <- vector("list", nrow(Xpost))
  logs <- vector("list", 2L * nrow(Xpost))
  counter <- 0L
  for (i in seq_len(nrow(Xpost))) {
    fiber <- trsc_evaluation_set(
      geometry, Xpost[i, ], representative = anchor
    )
    rows[[i]] <- data.frame(
      q = q_value, geometry_role = geometry_role, year = post$year[i],
      lower = fiber$interval[["lower"]], upper = fiber$interval[["upper"]],
      width = trsc_evaluation_width(fiber),
      point_identified = trsc_point_identified(fiber),
      status = fiber$status,
      lower_certificate = fiber$lower_result$certificate$method,
      upper_certificate = fiber$upper_result$certificate$method,
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
    logs[[counter]] <- support_record(
      fiber$lower_result, q_value, NA_real_, i, "fiber_lower", post$year[i]
    )
    counter <- counter + 1L
    logs[[counter]] <- support_record(
      fiber$upper_result, q_value, NA_real_, i, "fiber_upper", post$year[i]
    )
  }
  list(values = do.call(rbind, rows), logs = do.call(rbind, logs))
}

all_results <- list(
  run_fibers(primary_geometry, "primary_calibration", q1_anchor, 1L),
  run_fibers(primary_geometry, "primary_calibration", q2_boundary$level, 2L),
  run_fibers(full_pre_geometry, "full_pre_diagnostic", q1_anchor, 1L),
  run_fibers(full_pre_geometry, "full_pre_diagnostic", q2_boundary$level, 2L)
)
write.csv(do.call(rbind, lapply(all_results, `[[`, "values")),
          file.path(output_root, "tables", "evaluation_fibers.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(all_results, `[[`, "logs")),
          file.path(output_root, "solver-logs", "evaluation_fibers_solver.csv"), row.names = FALSE)
write_status("06_evaluation_fibers", "completed")
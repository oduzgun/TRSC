source(file.path("analysis", "california", "R", "00_helpers.R"))
library(TRSC)

q1 <- readRDS(file.path(output_root, "diagnostics", "q1_selection.rds"))
q2 <- readRDS(file.path(output_root, "diagnostics", "q2_selection.rds"))
q2_fit <- selected_fit(q2)
q2_h <- q2$selected$candidate_specification$h
boundary <- trsc_boundary_readouts(q2_fit, q2_h)

saveRDS(q1$selected, file.path(output_root, "diagnostics", "q1_boundary_level.rds"))
saveRDS(boundary, file.path(output_root, "diagnostics", "q2_boundary.rds"))
write.csv(data.frame(donor = names(boundary$level$representative),
                     level = boundary$level$representative,
                     previous = boundary$previous$representative),
          file.path(output_root, "diagnostics", "q2_boundary_representatives.csv"), row.names = FALSE)
write.csv(data.frame(coordinate = names(boundary$velocity_metric_coordinates),
                     velocity = boundary$velocity_metric_coordinates),
          file.path(output_root, "diagnostics", "q2_boundary_velocity.csv"), row.names = FALSE)
write.csv(data.frame(
  h = q2_h,
  level_cutoff = boundary$level$cutoff,
  previous_cutoff = boundary$previous$cutoff,
  velocity_norm = sqrt(sum(boundary$velocity_metric_coordinates^2)),
  stringsAsFactors = FALSE
), file.path(output_root, "tables", "q2_boundary_summary.csv"), row.names = FALSE)
write_status("05_boundary_readouts", "completed")
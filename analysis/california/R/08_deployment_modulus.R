source(file.path("analysis", "california", "R", "00_helpers.R"))
library(TRSC)

wide <- read.csv(processed_wide_path, check.names = FALSE)
donors <- setdiff(names(wide), c("year", "California"))
post <- wide[wide$year >= 1989, ]
Xpost <- as.matrix(post[donors])
colnames(Xpost) <- donors
primary <- readRDS(file.path(app_root, "data", "processed", "calibration_geometry.rds"))
full_pre <- trsc_geometry(as.matrix(wide[wide$year <= 1988, donors, drop = FALSE]))
epsilon_grid <- c(0, .25, .5, 1, 2, 5)
hull_width <- apply(Xpost, 1, max) - apply(Xpost, 1, min)

rows <- list(); logs <- list(); counter <- 0L
for (role in c("primary_calibration", "full_pre_diagnostic")) {
  geometry <- if (role == "primary_calibration") primary else full_pre
  for (i in seq_len(nrow(Xpost))) {
    for (epsilon in epsilon_grid) {
      ans <- trsc_deployment_modulus(geometry, Xpost[i, ], epsilon)
      counter <- counter + 1L
      rows[[counter]] <- data.frame(
        geometry_role = role, year = post$year[i], epsilon = epsilon,
        omega = ans$upper, hull_width = hull_width[i],
        hull_fraction = ans$upper / hull_width[i],
        exact_ambiguity = epsilon == 0 && ans$upper > 1e-7,
        status = ans$status, certificate_method = ans$certificate$method,
        stringsAsFactors = FALSE
      )
      logs[[counter]] <- support_record(ans, NA_integer_, NA_real_, i,
                                        paste0("deployment_epsilon_", epsilon), post$year[i])
    }
  }
}
write.csv(do.call(rbind, rows),
          file.path(output_root, "tables", "deployment_modulus.csv"), row.names = FALSE)
write.csv(do.call(rbind, logs),
          file.path(output_root, "solver-logs", "deployment_modulus_solver.csv"), row.names = FALSE)
write_status("08_deployment_modulus", "completed")
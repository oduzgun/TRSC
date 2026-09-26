source(file.path("analysis", "california", "R", "00_helpers.R"))
if (!requireNamespace("osqp", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) {
  stop("Static SCM requires the package dependencies osqp and Matrix.")
}

wide <- read.csv(processed_wide_path, check.names = FALSE)
donors <- setdiff(names(wide), c("year", "California"))

fit_static <- function(years, label) {
  rows <- match(years, wide$year)
  X <- as.matrix(wide[rows, donors, drop = FALSE])
  y <- wide$California[rows]
  P <- 2 * Matrix::Matrix(crossprod(X), sparse = TRUE)
  q <- -2 * drop(crossprod(X, y))
  A <- rbind(Matrix::Matrix(matrix(1, 1, length(donors)), sparse = TRUE),
             Matrix::Diagonal(length(donors)))
  result <- osqp::solve_osqp(
    P = P, q = q, A = A,
    l = c(1, rep(0, length(donors))),
    u = c(1, rep(1, length(donors))),
    pars = osqp::osqpSettings(
      verbose = FALSE, eps_abs = 1e-9, eps_rel = 1e-9,
      max_iter = 100000L, polishing = TRUE
    )
  )
  if (!tolower(trimws(result$info$status)) %in% c("solved", "solved inaccurate")) {
    stop(label, " static SCM failed: ", result$info$status)
  }
  w <- result$x
  names(w) <- donors
  pred <- drop(as.matrix(wide[donors]) %*% w)
  residual <- y - X %*% w
  write.csv(data.frame(donor = donors, weight = w),
            file.path(output_root, "tables", paste0(label, "_weights.csv")), row.names = FALSE)
  write.csv(data.frame(year = wide$year, observed = wide$California, synthetic = pred),
            file.path(output_root, "diagnostics", paste0(label, "_trajectory.csv")), row.names = FALSE)
  write.csv(data.frame(
    benchmark = label, start_year = min(years), end_year = max(years),
    rmse = sqrt(mean(residual^2)), mspe = mean(residual^2),
    solver_status = result$info$status, iterations = result$info$iter,
    primal_residual = result$info$prim_res, dual_residual = result$info$dual_res
  ), file.path(output_root, "solver-logs", paste0(label, "_solver.csv")), row.names = FALSE)
}

fit_static(1975:1985, "static_role_matched")
fit_static(1970:1988, "static_full_pre")

adh_status <- data.frame(
  benchmark = "original_adh_reproduction",
  status = if (requireNamespace("Synth", quietly = TRUE)) "available_not_automatically_run" else "missing_Synth_runtime",
  note = "Original predictor-based ADH reproduction is separate from outcome-only SCM.",
  stringsAsFactors = FALSE
)
write.csv(adh_status, file.path(output_root, "diagnostics", "adh_reproduction_status.csv"), row.names = FALSE)
write_status("03_static_scm", "completed")
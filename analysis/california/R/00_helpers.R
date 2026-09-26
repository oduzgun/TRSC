app_root <- file.path("analysis", "california")
output_root <- file.path(app_root, "output")

california_dirs <- c(
  file.path(app_root, "data", "raw"),
  file.path(app_root, "data", "processed"),
  file.path(output_root, "solver-logs"),
  file.path(output_root, "certificates"),
  file.path(output_root, "diagnostics"),
  file.path(output_root, "figures"),
  file.path(output_root, "tables")
)
invisible(lapply(california_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

processed_long_path <- file.path(app_root, "data", "processed", "smoking_long.csv")
processed_wide_path <- file.path(app_root, "data", "processed", "cigsale_wide.csv")

read_processed_long <- function() {
  x <- read.csv(processed_long_path, stringsAsFactors = FALSE, check.names = FALSE)
  x$year <- as.integer(x$year)
  x
}

selected_fit <- function(selection) {
  fit_name <- selection$selected$candidate_specification$fit_name
  fits <- attr(selection$candidates, "fits")
  if (is.null(fits) || is.null(fits[[fit_name]])) {
    stop("Selected frozen fit is not present in the candidate-library attributes.")
  }
  fits[[fit_name]]
}

support_record <- function(result, q, R, horizon, direction, year = NA_integer_) {
  data.frame(
    q = q,
    R = R,
    horizon = horizon,
    year = year,
    direction = direction,
    status = result$status,
    backend = result$backend,
    primal_objective = result$primal_objective,
    certified_upper = result$upper,
    primal_lower = result$primal_lower,
    certified_gap = result$certified_gap,
    primal_violation = result$primal_feasibility$violation,
    primal_exact = result$primal_feasibility$exact,
    certificate_method = result$certificate$method,
    stationarity_residual = result$certificate$stationarity_residual,
    stringsAsFactors = FALSE
  )
}

write_status <- function(stage, status, message = "") {
  write.csv(
    data.frame(stage = stage, status = status, message = message,
               stringsAsFactors = FALSE),
    file.path(output_root, "diagnostics", paste0(stage, "_status.csv")),
    row.names = FALSE
  )
}
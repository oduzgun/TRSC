if (!requireNamespace("TRSC", quietly = TRUE)) {
  stop(
    "Install the local package first from the repository root with: R CMD INSTALL .",
    call. = FALSE
  )
}

scripts <- sprintf(
  "%02d_%s.R",
  1:10,
  c(
    "retrieve_validate", "prepare_panel", "static_scm", "trsc_pre",
    "boundary_readouts", "evaluation_fibers", "dynamic_sets",
    "deployment_modulus", "outer_inference", "figures_tables"
  )
)

for (script in scripts) {
  message("Running ", script)
  source(file.path("analysis", "california", "R", script), echo = FALSE)
}

message("California workflow completed.")

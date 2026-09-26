source(file.path("analysis", "california", "R", "00_helpers.R"))

required <- c(
  "path_risk_envelope", "noise_covariance_bound", "restricted_coercivity_bound",
  "fitting_to_cutoff_transfer_bound", "predictive_to_causal_transfer_bound",
  "all_simplex_chord_epsilon", "valid_probability_allocation"
)
status <- data.frame(
  tier = c("tier_1_plugin", "tier_2_diagnostic_or_semi_feasible", "tier_3_honest_gate_5"),
  status = c("available_from_dynamic_sets", "not_run_without_supplied_bounds",
             "unavailable_missing_defensible_certificates"),
  confidence_interpretation = c(FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
write.csv(status, file.path(output_root, "certificates", "inference_tier_status.csv"), row.names = FALSE)
write.csv(data.frame(input = required, status = "not_supplied", stringsAsFactors = FALSE),
          file.path(output_root, "certificates", "missing_gate5_inputs.csv"), row.names = FALSE)
write_status(
  "09_outer_inference", "completed_without_honest_set",
  "No empirical rank, pseudo-oracle geometry discrepancy or tuning parameter was substituted for a Gate-5 certificate."
)
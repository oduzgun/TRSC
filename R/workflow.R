.trsc_validate_workflow_prepared <- function(prepared) {
  if (!inherits(prepared, "trsc_data") || !is.list(prepared) ||
      !identical(prepared$schema_version, 2L) ||
      !inherits(prepared$blocks, "trsc_blocks")) {
    stop("prepared must be a schema-version-2 trsc_data object.", call. = FALSE)
  }
  invisible(prepared)
}

.trsc_positive_integer <- function(x, label) {
  if (length(x) != 1L || !is.numeric(x) || is.complex(x) ||
      is.na(x) || !is.finite(x) || x != as.integer(x) || x < 1L) {
    stop(label, " must be one positive integer-valued number.", call. = FALSE)
  }
  as.integer(x)
}

.trsc_validate_sample_plan <- function(prepared, plan) {
  .trsc_validate_workflow_prepared(prepared)
  if (!inherits(plan, "trsc_sample_plan") || !is.list(plan) ||
      !identical(plan$schema_version, 1L) ||
      !identical(plan$schema, prepared$schema) ||
      !identical(plan$all_pre_periods, prepared$blocks$periods)) {
    stop("plan must be a compatible trsc_sample_plan for prepared.", call. = FALSE)
  }
  invisible(plan)
}

#' Construct the theorem-aligned pre-treatment sample-role plan
#'
#' @param prepared A schema-version-2 `trsc_data` object.
#' @param calibration_periods Contiguous raw pre-treatment calibration segment,
#'   strictly before the fitting block.
#' @param fitting_periods Contiguous intact fitting block of at least three
#'   dates.
#' @param validation_periods Contiguous raw validation segment after the fitting
#'   block and no later than the last untreated date.
#' @param calibration_thin Positive integer spacing for retained calibration
#'   observations. Retention begins at the first calibration date.
#' @param validation_thin Positive integer spacing for retained validation
#'   observations. Retention begins at the first validation date.
#' @param target_period Last untreated date. The current theorem-aligned
#'   implementation requires the final declared pre-treatment date.
#' @param min_gap0,min_gap1 Minimum separations, in declared period steps,
#'   between the last retained calibration date and the fitting block, and
#'   between the fitting block and first retained validation date.
#'
#' @return A `trsc_sample_plan` object containing raw and retained role
#'   assignments, effective gaps, counts, and deterministic thinning metadata.
#'
#' @details The fitting block is never thinned. Gaps are measured as date
#'   separations in units of `prepared$schema$period_step`, matching the
#'   dependence-transfer convention. The constructor verifies chronology and
#'   arithmetic only; it does not certify independence, beta-mixing, tail
#'   conditions, calibration complexity, or fitting assumptions.
#'
#' @export
trsc_sample_plan <- function(prepared,
                             calibration_periods,
                             fitting_periods,
                             validation_periods,
                             calibration_thin = 1L,
                             validation_thin = 1L,
                             target_period = NULL,
                             min_gap0 = 1L,
                             min_gap1 = 1L) {
  .trsc_validate_workflow_prepared(prepared)

  all_periods <- prepared$blocks$periods
  spacing <- prepared$schema$period_step
  if (is.null(target_period)) target_period <- max(all_periods)
  if (length(target_period) != 1L || !is.numeric(target_period) ||
      is.complex(target_period) || is.na(target_period) ||
      !is.finite(target_period) || target_period != max(all_periods)) {
    stop("target_period must equal the last declared pre-treatment date.",
         call. = FALSE)
  }

  calibration_thin <- .trsc_positive_integer(calibration_thin, "calibration_thin")
  validation_thin <- .trsc_positive_integer(validation_thin, "validation_thin")
  min_gap0 <- .trsc_positive_integer(min_gap0, "min_gap0")
  min_gap1 <- .trsc_positive_integer(min_gap1, "min_gap1")

  validate_segment <- function(x, label, minimum) {
    if (!is.numeric(x) || is.complex(x) || length(x) < minimum ||
        anyNA(x) || any(!is.finite(x)) || anyDuplicated(x) ||
        any(diff(x) <= 0)) {
      stop(label, " must contain at least ", minimum,
           " strictly increasing finite dates.", call. = FALSE)
    }
    index <- match(x, all_periods)
    if (anyNA(index)) {
      stop("Every ", label, " date must belong to the declared pre-treatment period.",
           call. = FALSE)
    }
    if (length(index) > 1L && any(diff(index) != 1L)) {
      stop(label, " must be one contiguous raw date segment.", call. = FALSE)
    }
    list(periods = x, index = as.integer(index))
  }

  calibration <- validate_segment(calibration_periods, "calibration_periods", 1L)
  fitting <- validate_segment(fitting_periods, "fitting_periods", 3L)
  validation <- validate_segment(validation_periods, "validation_periods", 2L)

  if (max(calibration$periods) >= min(fitting$periods)) {
    stop("The raw calibration segment must end before the fitting block.",
         call. = FALSE)
  }
  if (max(fitting$periods) >= min(validation$periods)) {
    stop("The fitting block must end before the validation segment.",
         call. = FALSE)
  }
  if (max(validation$periods) > target_period) {
    stop("The validation segment cannot extend beyond target_period.",
         call. = FALSE)
  }

  calibration_keep <- seq.int(1L, length(calibration$periods),
                              by = calibration_thin)
  validation_keep <- seq.int(1L, length(validation$periods),
                             by = validation_thin)
  calibration_retained <- calibration$periods[calibration_keep]
  validation_retained <- validation$periods[validation_keep]
  if (length(validation_retained) < 2L) {
    stop("validation_thin leaves fewer than two retained validation dates.",
         call. = FALSE)
  }

  gap0_raw <- (min(fitting$periods) - max(calibration_retained)) / spacing
  gap1_raw <- (min(validation_retained) - max(fitting$periods)) / spacing
  if (!is.finite(gap0_raw) || !is.finite(gap1_raw) ||
      gap0_raw != as.integer(gap0_raw) || gap1_raw != as.integer(gap1_raw)) {
    stop("Sample-role gaps are not aligned with the declared period_step.",
         call. = FALSE)
  }
  gap0 <- as.integer(gap0_raw)
  gap1 <- as.integer(gap1_raw)
  if (gap0 < min_gap0) {
    stop("The calibration-to-fitting separation is smaller than min_gap0.",
         call. = FALSE)
  }
  if (gap1 < min_gap1) {
    stop("The fitting-to-validation separation is smaller than min_gap1.",
         call. = FALSE)
  }

  expected_nc <- 1L + floor((length(calibration$periods) - 1L) /
                            calibration_thin)
  expected_nv <- 1L + floor((length(validation$periods) - 1L) /
                            validation_thin)

  structure(
    list(
      schema_version = 1L,
      schema = prepared$schema,
      all_pre_periods = all_periods,
      target_period = target_period,
      calibration = list(
        raw_periods = calibration$periods,
        raw_index = calibration$index,
        retained_periods = calibration_retained,
        retained_index = calibration$index[calibration_keep],
        raw_length = length(calibration$periods),
        retained_count = length(calibration_retained)
      ),
      fitting = list(
        periods = fitting$periods,
        index = fitting$index,
        length = length(fitting$periods),
        cutoff = max(fitting$periods),
        intact = TRUE
      ),
      validation = list(
        raw_periods = validation$periods,
        raw_index = validation$index,
        retained_periods = validation_retained,
        retained_index = validation$index[validation_keep],
        raw_length = length(validation$periods),
        retained_count = length(validation_retained)
      ),
      gaps = list(
        g0_steps = gap0,
        g1_steps = gap1,
        omitted_after_calibration = max(0L, gap0 - 1L),
        omitted_after_fitting = max(0L, gap1 - 1L)
      ),
      specification = list(
        calibration_thin = calibration_thin,
        validation_thin = validation_thin,
        min_gap0 = min_gap0,
        min_gap1 = min_gap1,
        deterministic = TRUE
      ),
      diagnostics = list(
        calibration_count_matches_formula =
          identical(as.integer(length(calibration_retained)),
                    as.integer(expected_nc)),
        validation_count_matches_formula =
          identical(as.integer(length(validation_retained)),
                    as.integer(expected_nv)),
        role_overlap = FALSE,
        fitting_block_intact = TRUE,
        theorem_assumptions_certified = FALSE
      )
    ),
    class = "trsc_sample_plan"
  )
}

#' Learn the donor prediction geometry from retained calibration dates
#'
#' @param prepared A schema-version-2 `trsc_data` object.
#' @param plan A compatible `trsc_sample_plan`.
#' @param tolerance Optional numerical rank tolerance passed to
#'   [trsc_geometry()].
#'
#' @return A `trsc_geometry` object with explicit calibration provenance.
#'
#' @details Only donor outcomes at the retained calibration dates are used.
#'   Target outcomes and validation outcomes are not used. The resulting
#'   empirical geometry remains a learned computational object and does not
#'   certify the population reference law.
#'
#' @export
trsc_calibration_geometry <- function(prepared, plan, tolerance = NULL) {
  .trsc_validate_sample_plan(prepared, plan)
  Xc <- prepared$blocks$Y0[plan$calibration$retained_index, , drop = FALSE]
  geometry <- trsc_geometry(Xc, tolerance = tolerance)
  geometry$provenance <- list(
    role = "calibration",
    raw_periods = plan$calibration$raw_periods,
    retained_periods = plan$calibration$retained_periods,
    thinning = plan$specification$calibration_thin,
    target_outcomes_used = FALSE,
    fitting_outcomes_used = FALSE,
    validation_outcomes_used = FALSE,
    population_reference_law_certified = FALSE,
    mixing_transfer_certified = FALSE
  )
  geometry
}

.trsc_validate_candidate_grid <- function(q, kappa) {
  if (!is.numeric(q) || is.complex(q) || !length(q) || anyNA(q) ||
      any(!is.finite(q)) || any(!(q %in% c(1, 2))) || anyDuplicated(q)) {
    stop("q must be a nonempty distinct vector containing only 1 and/or 2.",
         call. = FALSE)
  }
  if (!is.numeric(kappa) || is.complex(kappa) || !length(kappa) ||
      anyNA(kappa) || any(!is.finite(kappa)) || any(kappa <= 0) ||
      anyDuplicated(kappa)) {
    stop("kappa must contain distinct finite strictly positive values.",
         call. = FALSE)
  }
  list(q = as.integer(q), kappa = as.numeric(kappa))
}

#' Build the frozen finite TRSC fit/readout library
#'
#' @param prepared A schema-version-2 `trsc_data` object.
#' @param plan A compatible theorem-aligned `trsc_sample_plan`.
#' @param kappa Prespecified distinct positive regularization parameters.
#' @param q Prespecified difference orders, each equal to 1 or 2.
#' @param bandwidths Optional fixed terminal bandwidth library. If NULL, use the
#'   dyadic-plus-full-length library for each fit.
#' @param geometry Optional precomputed compatible `trsc_geometry`. If NULL,
#'   construct it from the retained calibration dates.
#' @param solver,solver_control,allow_inaccurate,feasibility_tolerance Passed to
#'   [trsc_fit()].
#'
#' @return A list-like `trsc_candidate_library` of frozen
#'   `trsc_readout` objects. Attributes retain the unique fitted paths,
#'   candidate grid, supplied geometry, sample plan, and deterministic tie rule.
#'
#' @details Each `(q, kappa)` fit is computed once on the intact fitting
#'   block. Readout bandwidths reuse that frozen fit and therefore do not create
#'   additional calibration events. Candidate order is fixed by the declared
#'   order of `q`, then `kappa`, then bandwidth. Exact numerical validation
#'   ties are resolved by the first candidate in this frozen order.
#'
#' @export
trsc_candidate_library <- function(prepared, plan, kappa,
                                   q = c(1L, 2L),
                                   bandwidths = NULL,
                                   geometry = NULL,
                                   solver = c("osqp"),
                                   solver_control = list(),
                                   allow_inaccurate = FALSE,
                                   feasibility_tolerance = 1e-6) {
  .trsc_validate_sample_plan(prepared, plan)
  grid_values <- .trsc_validate_candidate_grid(q, kappa)
  q <- grid_values$q
  kappa <- grid_values$kappa
  solver <- match.arg(solver)

  if (is.null(geometry)) {
    geometry <- trsc_calibration_geometry(prepared, plan)
  } else {
    geometry <- .trsc_align_geometry(
      geometry, colnames(prepared$blocks$Y0)
    )
    provenance <- geometry$provenance
    if (!is.list(provenance) ||
        !identical(provenance$role, "calibration") ||
        !identical(provenance$retained_periods,
                   plan$calibration$retained_periods) ||
        !identical(provenance$target_outcomes_used, FALSE) ||
        !identical(provenance$fitting_outcomes_used, FALSE) ||
        !identical(provenance$validation_outcomes_used, FALSE)) {
      stop(
        "A supplied geometry for trsc_candidate_library() must carry ",
        "calibration provenance matching plan. Use ",
        "trsc_calibration_geometry(prepared, plan) or leave geometry = NULL.",
        call. = FALSE
      )
    }
  }

  fits <- list()
  candidates <- list()
  grid_rows <- list()
  fit_counter <- 0L
  candidate_counter <- 0L

  for (q_value in q) {
    for (kappa_value in kappa) {
      fit_counter <- fit_counter + 1L
      fit_name <- sprintf("fit_%03d", fit_counter)
      fit <- trsc_fit(
        prepared = prepared,
        geometry = geometry,
        kappa = kappa_value,
        q = q_value,
        periods = plan$fitting$periods,
        solver = solver,
        solver_control = solver_control,
        allow_inaccurate = allow_inaccurate,
        feasibility_tolerance = feasibility_tolerance
      )
      fits[[fit_name]] <- fit
      readouts <- trsc_readout_library(fit, bandwidths = bandwidths)

      for (readout in readouts) {
        candidate_counter <- candidate_counter + 1L
        candidate_name <- sprintf("m%03d", candidate_counter)
        readout$candidate_specification <- list(
          candidate_index = candidate_counter,
          fit_name = fit_name,
          q = q_value,
          kappa = kappa_value,
          h = readout$bandwidth,
          tie_order = candidate_counter
        )
        candidates[[candidate_name]] <- readout
        grid_rows[[candidate_counter]] <- data.frame(
          candidate = candidate_name,
          fit = fit_name,
          q = q_value,
          kappa = kappa_value,
          h = readout$bandwidth,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  grid <- do.call(rbind, grid_rows)
  rownames(grid) <- NULL
  class(candidates) <- c("trsc_candidate_library", "list")
  attr(candidates, "schema_version") <- 1L
  attr(candidates, "grid") <- grid
  attr(candidates, "fits") <- fits
  attr(candidates, "geometry") <- geometry
  attr(candidates, "plan") <- plan
  attr(candidates, "frozen") <- TRUE
  attr(candidates, "tie_rule") <- "first exact minimum in frozen candidate order"
  candidates
}

#' Run the complete frozen pre-treatment TRSC workflow
#'
#' @param prepared A schema-version-2 `trsc_data` object.
#' @param plan A compatible `trsc_sample_plan`.
#' @param kappa,q,bandwidths Candidate-library specification passed to
#'   [trsc_candidate_library()].
#' @param geometry Optional precomputed calibration geometry.
#' @param kernel_weights Optional deterministic validation kernel weights passed
#'   to [trsc_select()].
#' @param mixing_beta Optional function returning an assumed absolute-regularity
#'   coefficient `beta(k)`. If supplied, the workflow reports the coupling
#'   bound from [trsc_mixing_bound()] but does not certify the assumption.
#' @param solver,solver_control,allow_inaccurate,feasibility_tolerance Passed to
#'   [trsc_candidate_library()].
#'
#' @return A `trsc_workflow` object containing the deterministic plan,
#'   calibration geometry, frozen candidate library, cutoff selection, selected
#'   readout, and optional mixing-bound diagnostic.
#'
#' @details Candidate specification is fully supplied before validation is
#'   scored. Calibration uses donor states only; fitting uses only the intact
#'   fitting block; validation is accessed only after all candidates are frozen.
#'   The function performs no rolling refit and no final full-sample refit.
#'   It is a deployable implementation of the paper's pre-treatment procedure,
#'   not a theorem-assumption verifier.
#'
#' @export
trsc_workflow <- function(prepared, plan, kappa,
                          q = c(1L, 2L),
                          bandwidths = NULL,
                          geometry = NULL,
                          kernel_weights = NULL,
                          mixing_beta = NULL,
                          solver = c("osqp"),
                          solver_control = list(),
                          allow_inaccurate = FALSE,
                          feasibility_tolerance = 1e-6) {
  .trsc_validate_sample_plan(prepared, plan)
  library <- trsc_candidate_library(
    prepared = prepared,
    plan = plan,
    kappa = kappa,
    q = q,
    bandwidths = bandwidths,
    geometry = geometry,
    solver = solver,
    solver_control = solver_control,
    allow_inaccurate = allow_inaccurate,
    feasibility_tolerance = feasibility_tolerance
  )
  learned_geometry <- attr(library, "geometry")

  selection <- trsc_select(
    prepared = prepared,
    candidates = library,
    validation_periods = plan$validation$retained_periods,
    target_period = plan$target_period,
    kernel_weights = kernel_weights
  )

  mixing <- NULL
  if (!is.null(mixing_beta)) {
    if (!is.function(mixing_beta)) {
      stop("mixing_beta must be NULL or a function beta(k).", call. = FALSE)
    }
    mixing <- trsc_mixing_bound(plan, mixing_beta)
  }

  structure(
    list(
      schema_version = 1L,
      call = match.call(),
      plan = plan,
      geometry = learned_geometry,
      candidate_library = library,
      selection = selection,
      selected = selection$selected,
      mixing = mixing,
      specification = list(
        q = as.integer(q),
        kappa = as.numeric(kappa),
        bandwidths = if (is.null(bandwidths)) "default_dyadic_plus_full"
                     else as.integer(bandwidths),
        candidate_order_fixed_before_validation = TRUE,
        tie_rule = attr(library, "tie_rule"),
        rolling_refit = FALSE,
        final_full_sample_refit = FALSE
      ),
      diagnostics = list(
        calibration_target_outcomes_used = FALSE,
        candidate_validation_outcomes_used = FALSE,
        validation_used_after_freeze = TRUE,
        theorem_assumptions_certified = FALSE,
        causal_identification = FALSE
      )
    ),
    class = "trsc_workflow"
  )
}

#' @export
print.trsc_sample_plan <- function(x, ...) {
  cat("TRSC pre-treatment sample plan\n")
  cat("  Calibration raw / retained: ", x$calibration$raw_length, " / ",
      x$calibration$retained_count, "\n", sep = "")
  cat("  Fitting dates (intact): ", x$fitting$length, "\n", sep = "")
  cat("  Validation raw / retained: ", x$validation$raw_length, " / ",
      x$validation$retained_count, "\n", sep = "")
  cat("  g0 / g1 separations: ", x$gaps$g0_steps, " / ",
      x$gaps$g1_steps, "\n", sep = "")
  cat("  Target period: ", x$target_period, "\n", sep = "")
  cat("  Theorem assumptions certified: no\n")
  invisible(x)
}

#' @export
print.trsc_candidate_library <- function(x, ...) {
  grid <- attr(x, "grid")
  cat("Frozen TRSC candidate library\n")
  cat("  Unique fitted paths: ", length(attr(x, "fits")), "\n", sep = "")
  cat("  Fit/readout candidates: ", length(x), "\n", sep = "")
  cat("  q values: ", paste(unique(grid$q), collapse = ", "), "\n", sep = "")
  cat("  Kappa values: ", paste(unique(grid$kappa), collapse = ", "), "\n", sep = "")
  cat("  Frozen before validation: yes\n")
  invisible(x)
}

#' @export
print.trsc_workflow <- function(x, ...) {
  cat("TRSC deployable pre-treatment workflow\n")
  cat("  Calibration retained: ", x$plan$calibration$retained_count, "\n", sep = "")
  cat("  Fitting dates: ", x$plan$fitting$length, "\n", sep = "")
  cat("  Validation retained: ", x$plan$validation$retained_count, "\n", sep = "")
  cat("  Candidates: ", length(x$candidate_library), "\n", sep = "")
  cat("  Selected: ", x$selection$selected_name, "\n", sep = "")
  cat("  Rolling refit: no\n")
  cat("  Final full-sample refit: no\n")
  cat("  Theorem assumptions certified: no\n")
  invisible(x)
}

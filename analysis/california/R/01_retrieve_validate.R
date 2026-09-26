source(file.path("analysis", "california", "R", "00_helpers.R"))

source_dta <- Sys.getenv("TRSC_CALIFORNIA_DTA", unset = "")
if (!nzchar(source_dta)) {
  candidates <- c(
    file.path(app_root, "data", "raw", "smoking.dta"),
    file.path("prism-uploads", "smoking.dta")
  )
  source_dta <- candidates[file.exists(candidates)][1L]
}

if (length(source_dta) && !is.na(source_dta) && file.exists(source_dta)) {
  if (requireNamespace("haven", quietly = TRUE)) {
    raw <- as.data.frame(haven::read_dta(source_dta))
    raw$state <- as.character(haven::as_factor(raw$state))
  } else if (requireNamespace("foreign", quietly = TRUE)) {
    raw <- foreign::read.dta(source_dta, convert.factors = TRUE)
    raw$state <- as.character(raw$state)
  } else {
    stop("Reading smoking.dta requires haven or the recommended foreign package.")
  }
} else if (requireNamespace("tidysynth", quietly = TRUE)) {
  data_env <- new.env(parent = emptyenv())
  utils::data("smoking", package = "tidysynth", envir = data_env)
  raw <- as.data.frame(data_env$smoking)
  raw$state <- as.character(raw$state)
} else {
  stop(
    "Supply smoking.dta via TRSC_CALIFORNIA_DTA or install the tidysynth package."
  )
}
raw$state <- trimws(raw$state)

required <- c("state", "year", "cigsale", "lnincome", "beer", "age15to24", "retprice")
if (!identical(names(raw), required)) stop("Unexpected smoking.dta schema.")
raw$year <- as.integer(raw$year)
if (nrow(raw) != 1209L) stop("Expected 1209 state-year rows.")
if (length(unique(raw$state)) != 39L) stop("Expected 39 states.")
if (!identical(range(raw$year), c(1970L, 2000L))) stop("Expected years 1970--2000.")
if (anyNA(raw$cigsale) || any(!is.finite(raw$cigsale))) stop("Cigarette sales must be complete.")
counts <- table(raw$state)
if (any(counts != 31L)) stop("Every state must have 31 annual observations.")
if (!"California" %in% raw$state) stop("California is absent.")

write.csv(raw, file.path(app_root, "data", "raw", "smoking.csv"), row.names = FALSE)
write.csv(
  data.frame(
    check = c("rows", "states", "year_min", "year_max", "missing_cigsale",
              "pre_periods", "post_periods", "donors"),
    value = c(nrow(raw), length(unique(raw$state)), min(raw$year), max(raw$year),
              sum(is.na(raw$cigsale)), sum(1970:1988 %in% raw$year),
              sum(1989:2000 %in% raw$year), 38),
    stringsAsFactors = FALSE
  ),
  file.path(output_root, "diagnostics", "data_validation.csv"),
  row.names = FALSE
)
write_status("01_retrieve_validate", "completed")

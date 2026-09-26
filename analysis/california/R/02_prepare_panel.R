source(file.path("analysis", "california", "R", "00_helpers.R"))

raw_path <- file.path(app_root, "data", "raw", "smoking.csv")
if (!file.exists(raw_path)) stop("Run 01_retrieve_validate.R first.")
raw <- read.csv(raw_path, stringsAsFactors = FALSE, check.names = FALSE)
raw$year <- as.integer(raw$year)

states <- sort(unique(raw$state))
donors <- setdiff(states, "California")
if (length(donors) != 38L) stop("The processed donor universe must contain 38 states.")

wide <- reshape(
  raw[c("year", "state", "cigsale")],
  idvar = "year", timevar = "state", direction = "wide"
)
names(wide) <- sub("^cigsale\\.", "", names(wide))
wide <- wide[order(wide$year), c("year", "California", donors)]
if (anyNA(wide)) stop("Processed cigarette-sales panel is not rectangular.")

write.csv(raw, processed_long_path, row.names = FALSE)
write.csv(wide, processed_wide_path, row.names = FALSE)
write.csv(data.frame(donor = donors),
          file.path(app_root, "data", "processed", "donors.csv"), row.names = FALSE)
write_status("02_prepare_panel", "completed")
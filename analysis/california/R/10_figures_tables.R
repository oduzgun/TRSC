source(file.path("analysis", "california", "R", "00_helpers.R"))

wide <- read.csv(processed_wide_path, check.names = FALSE)
role <- read.csv(file.path(output_root, "diagnostics", "static_role_matched_trajectory.csv"))
full <- read.csv(file.path(output_root, "diagnostics", "static_full_pre_trajectory.csv"))
q1_path <- read.csv(file.path(output_root, "diagnostics", "q1_fitted_path.csv"))
q2_path <- read.csv(file.path(output_root, "diagnostics", "q2_fitted_path.csv"))

png(file.path(output_root, "figures", "figure1_static_and_observed.png"), 1800, 1100, res = 180)
plot(wide$year, wide$California, type = "l", lwd = 2, col = "black",
     xlab = "Year", ylab = "Cigarette sales", main = "California and static SCM benchmarks")
lines(role$year, role$synthetic, col = "#2C7FB8", lwd = 2)
lines(full$year, full$synthetic, col = "#D95F0E", lwd = 2, lty = 2)
abline(v = 1989, lty = 3)
legend("topright", c("California observed", "Role-matched SCM", "Full-pre SCM"),
       col = c("black", "#2C7FB8", "#D95F0E"), lty = c(1, 1, 2), lwd = 2, bty = "n")
dev.off()

png(file.path(output_root, "figures", "figure2_trsc_pre_paths.png"), 1800, 1000, res = 180)
yr <- sort(unique(c(q1_path$period, q2_path$period)))
ylim <- range(c(q1_path$observed, q1_path$fitted, q2_path$fitted))
plot(q1_path$period, q1_path$observed, type = "l", lwd = 2, col = "black", ylim = ylim,
     xlab = "Fitting year", ylab = "Cigarette sales", main = "Selected TRSC fitting paths")
lines(q1_path$period, q1_path$fitted, col = "#1B9E77", lwd = 2)
lines(q2_path$period, q2_path$fitted, col = "#7570B3", lwd = 2, lty = 2)
legend("topright", c("California", "q=1 fitted", "q=2 fitted"),
       col = c("black", "#1B9E77", "#7570B3"), lty = c(1, 1, 2), lwd = 2, bty = "n")
dev.off()

plot_widths <- function(q_value) {
  x <- read.csv(file.path(output_root, "tables", paste0("q", q_value, "_dynamic_bounds.csv")))
  png(file.path(output_root, "figures", paste0("figure_q", q_value, "_width_vs_R.png")), 1800, 1100, res = 180)
  matplot(sort(unique(x$R)), sapply(c(1,3,6,12), function(h) x$width[x$horizon == h]),
          type = "l", lty = 1, lwd = 2, xlab = "R", ylab = "Identified-set width",
          main = paste0("q=", q_value, " width sensitivity"))
  legend("topleft", paste("h=", c(1,3,6,12), sep = ""), col = 1:4, lty = 1, lwd = 2, bty = "n")
  dev.off()
}
plot_widths(1L); plot_widths(2L)

write_status("10_figures_tables", "completed")
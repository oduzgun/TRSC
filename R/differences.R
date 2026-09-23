# Internal temporal-difference infrastructure for TRSC.

.trsc_difference_matrix <- function(n, q) {
  if (length(n) != 1L || !is.numeric(n) || is.complex(n) || is.na(n) ||
      !is.finite(n) || n != as.integer(n) || n < 2L) {
    stop("n must be one integer-valued number of at least 2.", call. = FALSE)
  }
  if (length(q) != 1L || !is.numeric(q) || is.complex(q) || is.na(q) ||
      !is.finite(q) || !(q %in% c(1, 2))) {
    stop("q must be exactly 1 or 2.", call. = FALSE)
  }
  n <- as.integer(n)
  q <- as.integer(q)
  if (n <= q) stop("n must exceed q.", call. = FALSE)

  D <- matrix(0, nrow = n - q, ncol = n)
  if (q == 1L) {
    for (i in seq_len(n - 1L)) D[i, i:(i + 1L)] <- c(-1, 1)
  } else {
    for (i in seq_len(n - 2L)) D[i, i:(i + 2L)] <- c(1, -2, 1)
  }
  D
}

.trsc_temporal_energy <- function(weights, geometry, q) {
  if (!inherits(geometry, "trsc_geometry")) {
    stop("geometry must be a trsc_geometry object.", call. = FALSE)
  }
  if (!is.matrix(weights) || !is.numeric(weights) || is.complex(weights) ||
      any(!is.finite(weights)) || nrow(weights) != geometry$donors) {
    stop("weights must be a finite real donor-by-time numeric matrix.", call. = FALSE)
  }
  D <- .trsc_difference_matrix(ncol(weights), q)
  delta <- D %*% t(weights)
  sum((delta %*% geometry$Gamma) * delta)
}

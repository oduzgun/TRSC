#' Construct empirical TRSC predictive geometry
#'
#' Construct the uncentered empirical second-moment geometry of donor states.
#' The matrix is a prediction metric, not a covariance matrix. Rank deficiency
#' is permitted and is reported rather than treated as an error.
#'
#' @param donor_states A finite numeric matrix or numeric data frame with
#'   observations in rows and donors in columns. At least two donor columns are
#'   required.
#' @param tolerance Optional nonnegative absolute eigenvalue tolerance used only
#'   for numerical rank diagnostics and construction of quotient coordinates.
#'   If NULL, a scale-aware machine-precision tolerance is used.
#'
#' @return An object of class `trsc_geometry` containing the empirical
#'   second-moment matrix `Gamma`, a factor `B` satisfying
#'   `crossprod(B) = Gamma` up to numerical precision on retained directions,
#'   and empirical rank diagnostics.
#'
#' @details No inverse of `Gamma` is formed. `representation_unique_empirical`
#'   is only an empirical finite-resolution diagnostic: it asks whether the
#'   retained quotient map is injective on the affine simplex. It is not a
#'   population identification certificate.
#'
#' @export
trsc_geometry <- function(donor_states, tolerance = NULL) {
  if (is.data.frame(donor_states)) {
    numeric_cols <- vapply(donor_states, function(x) is.numeric(x) && !is.complex(x), logical(1))
    if (!all(numeric_cols)) {
      stop("donor_states data-frame columns must all be real numeric.", call. = FALSE)
    }
    donor_states <- as.matrix(donor_states)
  }
  if (!is.matrix(donor_states) || !is.numeric(donor_states) || is.complex(donor_states)) {
    stop("donor_states must be a real numeric matrix or numeric data frame.", call. = FALSE)
  }
  if (nrow(donor_states) < 1L || ncol(donor_states) < 2L) {
    stop("donor_states must contain at least one observation and two donors.", call. = FALSE)
  }
  if (any(!is.finite(donor_states))) {
    stop("donor_states must contain only finite values.", call. = FALSE)
  }

  X <- unclass(donor_states)
  storage.mode(X) <- "double"
  J <- ncol(X)
  donor_names <- colnames(X)
  if (is.null(donor_names)) donor_names <- paste0("donor_", seq_len(J))
  if (anyNA(donor_names) || any(!nzchar(donor_names)) || anyDuplicated(donor_names)) {
    stop("donor column names, when supplied, must be unique and nonempty.", call. = FALSE)
  }
  colnames(X) <- donor_names

  Gamma <- crossprod(X) / nrow(X)
  Gamma <- (Gamma + t(Gamma)) / 2
  dimnames(Gamma) <- list(donor_names, donor_names)

  eig <- eigen(Gamma, symmetric = TRUE)
  spectral_scale <- max(c(1, abs(eig$values)))
  if (is.null(tolerance)) {
    tolerance <- max(dim(X)) * .Machine$double.eps * spectral_scale
  } else {
    if (length(tolerance) != 1L || !is.numeric(tolerance) || is.complex(tolerance) ||
        is.na(tolerance) || !is.finite(tolerance) || tolerance < 0) {
      stop("tolerance must be NULL or one finite nonnegative real number.", call. = FALSE)
    }
  }
  if (any(eig$values < -100 * tolerance - 100 * .Machine$double.eps * spectral_scale)) {
    stop("Numerical failure: empirical second-moment matrix is not positive semidefinite.", call. = FALSE)
  }

  keep <- eig$values > tolerance
  rank <- sum(keep)
  if (rank > 0L) {
    vectors <- eig$vectors[, keep, drop = FALSE]
    values <- eig$values[keep]
    B <- diag(sqrt(values), nrow = rank) %*% t(vectors)
  } else {
    B <- matrix(numeric(0), nrow = 0L, ncol = J)
  }
  colnames(B) <- donor_names
  rownames(B) <- if (rank) paste0("profile_", seq_len(rank)) else character(0)

  affine_map <- rbind(B, rep(1, J))
  affine_rank <- qr(affine_map, tol = sqrt(.Machine$double.eps))$rank
  representation_unique <- identical(as.integer(affine_rank), as.integer(J))

  structure(
    list(
      schema_version = 1L,
      Gamma = Gamma,
      B = B,
      donor_names = donor_names,
      observations = nrow(X),
      donors = J,
      tolerance = tolerance,
      diagnostics = list(
        empirical_rank = as.integer(rank),
        empirical_null_dimension = as.integer(J - rank),
        affine_simplex_rank = as.integer(affine_rank),
        representation_unique_empirical = representation_unique,
        eigenvalues = eig$values,
        centered = FALSE
      )
    ),
    class = "trsc_geometry"
  )
}

#' Empirical distance between two donor representations
#'
#' @param geometry A `trsc_geometry` object.
#' @param w,v Finite numeric donor-coordinate vectors aligned with
#'   `geometry$donor_names`.
#'
#' @return The empirical predictive distance induced by `geometry$Gamma`.
#' @export
trsc_profile_distance <- function(geometry, w, v) {
  if (!inherits(geometry, "trsc_geometry") || !is.list(geometry) ||
      is.null(geometry$Gamma) || is.null(geometry$donor_names)) {
    stop("geometry must be a trsc_geometry object.", call. = FALSE)
  }
  align <- function(x, label) {
    if (!is.numeric(x) || is.complex(x) || length(x) != geometry$donors || any(!is.finite(x))) {
      stop(label, " must be a finite real numeric vector with one value per donor.", call. = FALSE)
    }
    if (!is.null(names(x))) {
      if (anyNA(names(x)) || anyDuplicated(names(x)) ||
          !setequal(names(x), geometry$donor_names)) {
        stop(label, " names must cover geometry donor names exactly.", call. = FALSE)
      }
      x <- x[geometry$donor_names]
    }
    as.numeric(x)
  }
  w <- align(w, "w")
  v <- align(v, "v")
  d <- w - v
  sq <- drop(crossprod(d, geometry$Gamma %*% d))
  numerical_scale <- max(1, sum(abs(d))^2 * max(abs(geometry$Gamma)))
  if (sq < -100 * .Machine$double.eps * numerical_scale) {
    stop("Numerical failure while evaluating predictive distance.", call. = FALSE)
  }
  sqrt(max(0, sq))
}

#' Test empirical predictive-profile equivalence
#'
#' @param geometry A `trsc_geometry` object.
#' @param w,v Donor-coordinate vectors.
#' @param tolerance Finite nonnegative distance tolerance.
#' @return A single logical value.
#' @export
trsc_profile_equivalent <- function(geometry, w, v, tolerance = sqrt(.Machine$double.eps)) {
  if (length(tolerance) != 1L || !is.numeric(tolerance) || is.complex(tolerance) ||
      is.na(tolerance) || !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be one finite nonnegative real number.", call. = FALSE)
  }
  trsc_profile_distance(geometry, w, v) <= tolerance
}

#' Quotient coordinates of a donor representation
#'
#' @param geometry A `trsc_geometry` object.
#' @param w A finite donor-coordinate vector.
#' @return Empirical quotient coordinates `B %*% w`. These coordinates
#'   describe the fitted-resolution predictive profile, not structural donor
#'   shares and not a causal counterfactual.
#' @export
trsc_profile_coordinates <- function(geometry, w) {
  if (!inherits(geometry, "trsc_geometry")) {
    stop("geometry must be a trsc_geometry object.", call. = FALSE)
  }
  if (!is.numeric(w) || is.complex(w) || length(w) != geometry$donors || any(!is.finite(w))) {
    stop("w must be a finite real numeric vector with one value per donor.", call. = FALSE)
  }
  if (!is.null(names(w))) {
    if (anyNA(names(w)) || anyDuplicated(names(w)) ||
        !setequal(names(w), geometry$donor_names)) {
      stop("w names must cover geometry donor names exactly.", call. = FALSE)
    }
    w <- w[geometry$donor_names]
  }
  drop(geometry$B %*% as.numeric(w))
}

#' @export
print.trsc_geometry <- function(x, ...) {
  cat("TRSC empirical predictive geometry\n")
  cat("  Observations: ", x$observations, "\n", sep = "")
  cat("  Donors:       ", x$donors, "\n", sep = "")
  cat("  Empirical rank: ", x$diagnostics$empirical_rank, "\n", sep = "")
  cat("  Null dimension: ", x$diagnostics$empirical_null_dimension, "\n", sep = "")
  cat("  Raw representation unique at empirical resolution: ",
      if (isTRUE(x$diagnostics$representation_unique_empirical)) "yes" else "no", "\n", sep = "")
  invisible(x)
}

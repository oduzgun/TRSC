# Internal convex-optimization infrastructure for TRSC.

.trsc_align_geometry <- function(geometry, donor_names) {
  if (!inherits(geometry, "trsc_geometry") || !is.list(geometry) ||
      is.null(geometry$Gamma) || is.null(geometry$B) ||
      is.null(geometry$donor_names)) {
    stop("geometry must be a trsc_geometry object.", call. = FALSE)
  }
  donor_names <- as.character(donor_names)
  if (length(donor_names) != geometry$donors ||
      !setequal(donor_names, geometry$donor_names)) {
    stop("geometry donor names must match the fitted donor set exactly.", call. = FALSE)
  }
  if (identical(donor_names, geometry$donor_names)) return(geometry)

  order <- match(donor_names, geometry$donor_names)
  out <- geometry
  out$Gamma <- geometry$Gamma[order, order, drop = FALSE]
  out$B <- geometry$B[, order, drop = FALSE]
  out$donor_names <- donor_names
  dimnames(out$Gamma) <- list(donor_names, donor_names)
  colnames(out$B) <- donor_names
  out
}

.trsc_build_qp <- function(y, X, geometry, kappa, q) {
  n <- length(y)
  J <- ncol(X)
  D <- .trsc_difference_matrix(n, q)

  fit_blocks <- lapply(seq_len(n), function(i) {
    Matrix::Matrix(tcrossprod(X[i, ]), sparse = TRUE)
  })
  H_fit <- do.call(Matrix::bdiag, fit_blocks)

  DtD <- Matrix::Matrix(crossprod(D), sparse = TRUE)
  Gamma <- Matrix::Matrix(geometry$Gamma, sparse = TRUE)
  H_penalty <- kronecker(DtD, Gamma)

  P <- 2 * (H_fit + kappa * H_penalty)
  P <- Matrix::forceSymmetric(P, uplo = "U")
  qvec <- -2 * as.vector(t(X * y))

  sum_row <- Matrix::Matrix(matrix(1, nrow = 1L, ncol = J), sparse = TRUE)
  simplex_equalities <- kronecker(Matrix::Diagonal(n), sum_row)
  coordinate_bounds <- Matrix::Diagonal(n * J)
  A <- rbind(simplex_equalities, coordinate_bounds)

  l <- c(rep(1, n), rep(0, n * J))
  u <- c(rep(1, n), rep(1, n * J))

  list(P = P, q = qvec, A = A, l = l, u = u, D = D)
}

.trsc_solve_osqp <- function(problem, solver_control = list(),
                             allow_inaccurate = FALSE) {
  if (!is.list(solver_control) || is.null(names(solver_control)) &&
      length(solver_control)) {
    stop("solver_control must be a named list.", call. = FALSE)
  }

  defaults <- list(
    verbose = FALSE,
    eps_abs = 1e-8,
    eps_rel = 1e-8,
    max_iter = 100000L,
    polishing = TRUE
  )
  settings_args <- utils::modifyList(defaults, solver_control)
  settings <- do.call(osqp::osqpSettings, settings_args)

  result <- osqp::solve_osqp(
    P = problem$P,
    q = problem$q,
    A = problem$A,
    l = problem$l,
    u = problem$u,
    pars = settings
  )

  status <- tolower(trimws(as.character(result$info$status)))
  solved <- identical(status, "solved")
  solved_inaccurate <- identical(status, "solved inaccurate")
  if (!solved && !(allow_inaccurate && solved_inaccurate)) {
    stop(
      "TRSC optimization did not return an accepted solution (OSQP status: ",
      result$info$status, ").",
      call. = FALSE
    )
  }
  if (is.null(result$x) || any(!is.finite(result$x))) {
    stop("TRSC optimization returned a nonfinite primal solution.", call. = FALSE)
  }

  result
}

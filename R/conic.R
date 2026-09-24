# Conic support optimization and exact rational weak-duality certificates.
# The mathematical affine coefficients are retained as bigq rationals. ECOS
# receives double approximations, but its status/point estimate is not a proof.

.trsc_Q <- function(x) if (inherits(x, "bigq")) x else gmp::as.bigq(x)
.trsc_mm <- function(a, b) gmp::`%*%`(.trsc_Q(a), .trsc_Q(b))
.trsc_num <- function(x) if (inherits(x, "bigq")) gmp::asNumeric(x) else x

.trsc_up <- function(x) {
  x <- .trsc_Q(x)
  v <- as.numeric(.trsc_num(x))
  if (length(v) != 1L || !is.finite(v)) stop("Nonfinite numerical certificate.", call. = FALSE)
  while (.trsc_Q(v) < x) {
    v <- v + max(abs(v) * 4 * .Machine$double.eps, .Machine$double.xmin)
    if (!is.finite(v)) stop("Certificate exceeds finite numeric range.", call. = FALSE)
  }
  v
}

.trsc_norm_up <- function(x) {
  x <- .trsc_Q(x)
  if (!length(x) || all(x == 0)) return(0)
  scale <- .trsc_up(max(abs(x)))
  s <- sum(x * x)
  v <- scale * sqrt(sum(as.numeric(.trsc_num(x / .trsc_Q(scale)))^2))
  if (!is.finite(v)) stop("Norm exceeds finite numeric range.", call. = FALSE)
  while (.trsc_Q(v)^2 < s) v <- v + max(abs(v) * 4 * .Machine$double.eps, .Machine$double.xmin)
  v
}

.trsc_scalar <- function(x, label, lower = 0, upper = Inf) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L ||
      !is.finite(x) || x < lower || x > upper)
    stop(label, " must be one finite number in [", lower, ", ", upper, "].", call. = FALSE)
  as.numeric(x)
}

.trsc_conic_problem <- function(n, blocks, extra_upper = numeric()) {
  n <- as.integer(n)
  covered <- unlist(blocks, use.names = FALSE)
  if (anyDuplicated(covered) || any(covered < 1 | covered > n)) stop("Invalid simplex blocks.")
  extra <- setdiff(seq_len(n), covered)
  if (length(extra) != length(extra_upper)) stop("Extra variables require finite bounds.")
  A <- matrix(0, length(blocks), n)
  for (i in seq_along(blocks)) A[i, blocks[[i]]] <- 1
  list(n = n, blocks = blocks, extra = extra, extra_upper = .trsc_Q(extra_upper),
       A = .trsc_Q(A), b = .trsc_Q(rep(1, length(blocks))),
       linear = list(list(G = .trsc_Q(-diag(n)), h = .trsc_Q(rep(0, n)))),
       soc = list())
}

# Equalities apart from simplex sums use two inequalities. This avoids rank
# assumptions on profile/auxiliary equations and does not discard any equation.
.trsc_add_equality <- function(p, M, rhs) {
  M <- .trsc_Q(M); rhs <- .trsc_Q(rhs)
  if (ncol(M) != p$n || nrow(M) != length(rhs)) stop("Equality dimensions differ.")
  p$linear <- c(p$linear, list(list(G = M, h = rhs), list(G = -M, h = -rhs)))
  p
}

.trsc_add_soc <- function(p, M, offset, radius) {
  M <- .trsc_Q(M); offset <- .trsc_Q(offset); radius <- .trsc_Q(radius)
  if (ncol(M) != p$n || nrow(M) != length(offset) || length(radius) != 1 || radius < 0)
    stop("Invalid SOC dimensions or radius.")
  p$soc <- c(p$soc, list(list(G = rbind(.trsc_Q(matrix(0, 1, p$n)), -M),
                             h = c(radius, offset))))
  p
}

.trsc_pack_conic <- function(p) {
  linear <- p$linear
  if (length(p$extra)) {
    U <- matrix(0, length(p$extra), p$n)
    U[cbind(seq_along(p$extra), p$extra)] <- 1
    linear <- c(linear, list(list(G = .trsc_Q(U), h = p$extra_upper)))
  }
  parts <- c(linear, p$soc)
  p$G <- do.call(rbind, lapply(parts, `[[`, "G"))
  p$h <- do.call(c, lapply(parts, `[[`, "h"))
  p$dims <- list(l = as.integer(sum(vapply(linear, function(v) nrow(v$G), integer(1)))),
                 q = as.integer(vapply(p$soc, function(v) length(v$h), integer(1))), e = 0L)
  p
}

.trsc_residual_support <- function(p, residual) {
  ans <- .trsc_Q(0)
  for (b in p$blocks) ans <- ans + max(residual[b])
  for (i in seq_along(p$extra)) ans <- ans + max(.trsc_Q(0), residual[p$extra[i]]) * p$extra_upper[i]
  ans
}

.trsc_certify_upper <- function(p, objective, y = NULL, z = NULL) {
  if (is.null(p$G)) p <- .trsc_pack_conic(p)
  objective <- .trsc_Q(objective)
  fallback <- .trsc_residual_support(p, objective)
  method <- "simplex_hull"
  bound <- fallback; residual_size <- NA_real_; dual_value <- NA_real_
  if (length(y) == nrow(p$A) && length(z) == nrow(p$G) &&
      all(is.finite(y)) && all(is.finite(z))) {
    z <- .trsc_Q(z); y <- .trsc_Q(y)
    if (p$dims$l) {
      negative <- seq_len(p$dims$l)[z[seq_len(p$dims$l)] < 0]
      if (length(negative)) z[negative] <- .trsc_Q(0)
    }
    first <- p$dims$l + 1L
    for (size in p$dims$q) {
      ind <- first + seq_len(size) - 1L
      z[ind[1]] <- max(z[ind[1]], .trsc_Q(.trsc_norm_up(z[ind[-1]])))
      stopifnot(z[ind[1]] >= 0, z[ind[1]]^2 >= sum(z[ind[-1]]^2))
      first <- first + size
    }
    residual <- objective - .trsc_mm(t(p$A), y) - .trsc_mm(t(p$G), z)
    dual <- sum(p$b * y) + sum(p$h * z)
    corrected <- dual + .trsc_residual_support(p, residual)
    if (corrected <= bound) { bound <- corrected; method <- "rational_dual_residual" }
    residual_size <- .trsc_up(max(abs(residual)))
    dual_value <- .trsc_up(dual)
  }
  list(valid = TRUE, upper = .trsc_up(bound), method = method,
       arithmetic = "exact rational affine arithmetic; checked SOC dual; outward conversion",
       stationarity_residual = residual_size, dual_objective = dual_value,
       scope = "supplied mathematical program; statistical assumptions are separate")
}

.trsc_primal_diagnostics <- function(p, x) {
  if (length(x) != p$n || any(!is.finite(x))) return(list(violation = Inf, exact = FALSE))
  xq <- .trsc_Q(x)
  slack <- p$h - .trsc_mm(p$G, xq)
  violation <- max(abs(as.numeric(.trsc_num(.trsc_mm(p$A, xq) - p$b))))
  exact <- all(.trsc_mm(p$A, xq) == p$b)
  if (p$dims$l) {
    lin <- slack[seq_len(p$dims$l)]
    violation <- max(violation, -.trsc_num(min(lin)))
    exact <- exact && all(lin >= 0)
  }
  first <- p$dims$l + 1L
  for (size in p$dims$q) {
    v <- slack[first + seq_len(size) - 1L]
    violation <- max(violation, .trsc_norm_up(v[-1]) - .trsc_num(v[1]))
    exact <- exact && v[1] >= 0 && v[1]^2 >= sum(v[-1]^2)
    first <- first + size
  }
  list(violation = max(0, violation), exact = isTRUE(exact))
}

.trsc_solve_socp <- function(p, objective, control = list()) {
  p <- .trsc_pack_conic(p)
  cq <- .trsc_Q(objective)
  if (length(cq) != p$n) stop("Objective dimension differs.")
  defaults <- list(maxit = 200L, feastol = 1e-9, reltol = 1e-9, abstol = 1e-9, verbose = 0L)
  if (!is.list(control) || (length(control) && (is.null(names(control)) || anyDuplicated(names(control)))))
    stop("control must be a uniquely named list.", call. = FALSE)
  settings <- do.call(ECOSolveR::ecos.control, utils::modifyList(defaults, control))
  raw <- tryCatch(ECOSolveR::ECOS_csolve(
    c = -as.numeric(.trsc_num(cq)),
    G = Matrix::Matrix(.trsc_num(p$G), sparse = TRUE), h = as.numeric(.trsc_num(p$h)),
    dims = p$dims, A = Matrix::Matrix(.trsc_num(p$A), sparse = TRUE),
    b = as.numeric(.trsc_num(p$b)), control = settings), error = function(e) list(error = conditionMessage(e)))
  cert <- .trsc_certify_upper(p, cq, raw$y, raw$z)
  infeas <- .trsc_certify_upper(p, rep(0, p$n), raw$y, raw$z)
  infeasible <- infeas$upper < 0
  d <- .trsc_primal_diagnostics(p, raw$x)
  value <- if (length(raw$x) == p$n && all(is.finite(raw$x))) sum(as.numeric(.trsc_num(cq)) * raw$x) else NA_real_
  lower <- if (d$exact) -.trsc_up(-sum(cq * .trsc_Q(raw$x))) else NA_real_
  list(backend = "ECOSolveR", status = if (infeasible) "structural_infeasible" else if (is.null(raw$error)) as.character(raw$infostring) else "solver_error",
       upper = if (infeasible) -Inf else cert$upper, primal_objective = value,
       primal_lower = lower, dual_objective = cert$dual_objective,
       primal_feasibility = d, dual_feasibility = cert$valid,
       certified_gap = if (is.finite(lower)) cert$upper - lower else NA_real_,
       certificate = cert, infeasibility_certificate = if (infeasible) infeas else NULL,
       tolerance = settings, solution = raw$x, raw = raw, program = p)
}

.trsc_solve_lp <- function(p, objective, control = list()) {
  if (length(p$soc)) stop("LP contains SOC constraints.")
  .trsc_solve_socp(p, objective, control)
}

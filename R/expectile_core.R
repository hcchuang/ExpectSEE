check_tau <- function(tau) {
  if (!is.numeric(tau) || length(tau) != 1 || tau <= 0 || tau >= 1) {
    stop("tau must be one number in (0, 1)")
  }
  invisible(tau)
}

safe_inv <- function(a, tol = 1e-10) {
  out <- tryCatch(solve(a), error = function(e) NULL)
  if (!is.null(out)) return(out)
  eig <- eigen((a + t(a)) / 2, symmetric = TRUE)
  eig$vectors %*% diag(ifelse(eig$values > tol, 1 / eig$values, 0)) %*% t(eig$vectors)
}

check_fixed_design <- function(x, label = "fixed-dimensional model") {
  if (ncol(x) >= nrow(x)) {
    stop(label, " is not estimable: p must be smaller than n (p = ",
         ncol(x), ", n = ", nrow(x), ")", call. = FALSE)
  }
  invisible(TRUE)
}

sym <- function(a) (a + t(a)) / 2

expectile_weight <- function(e, tau) {
  check_tau(tau)
  abs(tau - as.numeric(e < 0))
}

rho_dot <- function(e, tau) 2 * expectile_weight(e, tau) * e

omega_tau <- function(e, tau) 2 * expectile_weight(e, tau)

expectile_loss <- function(resid, tau) {
  expectile_weight(resid, tau) * resid^2
}

fit_expectile <- function(y, x, tau, beta0 = NULL, max_iter = 100, tol = 1e-13) {
  check_tau(tau)
  x <- as.matrix(x)
  y <- as.numeric(y)
  if (nrow(x) != length(y)) stop("x and y have incompatible dimensions")

  if (is.null(beta0)) {
    beta0 <- as.numeric(qr.coef(qr(x), y))
    if (any(!is.finite(beta0))) beta0 <- as.numeric(safe_inv(crossprod(x)) %*% crossprod(x, y))
  }

  beta <- as.numeric(beta0)
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    resid <- as.numeric(y - x %*% beta)
    w <- expectile_weight(resid, tau)
    a <- crossprod(x, x * w)
    b <- crossprod(x, y * w)
    beta_new <- as.numeric(tryCatch(solve(a, b), error = function(e) safe_inv(a) %*% b))
    if (max(abs(beta_new - beta)) < tol) {
      beta <- beta_new
      converged <- TRUE
      break
    }
    beta <- beta_new
  }

  fitted <- as.numeric(x %*% beta)
  resid <- y - fitted
  score <- colMeans(x * as.numeric(rho_dot(resid, tau)))
  structure(
    list(
      beta = beta, fitted = fitted, resid = resid,
      weights = expectile_weight(resid, tau),
      loss = expectile_loss(resid, tau),
      iter = iter, converged = converged,
      score_max = max(abs(score)), tau = tau,
      n = length(y), p = ncol(x)
    ),
    class = "expectsee_als"
  )
}

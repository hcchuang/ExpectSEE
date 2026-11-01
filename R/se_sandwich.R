meat_iid <- function(l_mat) crossprod(l_mat) / nrow(l_mat)

meat_hac <- function(l_mat, lag = 0) {
  n <- nrow(l_mat)
  omega <- crossprod(l_mat) / n
  if (lag <= 0) return(omega)
  for (h in seq_len(lag)) {
    g <- crossprod(l_mat[(h + 1):n, , drop = FALSE],
                   l_mat[1:(n - h), , drop = FALSE]) / n
    omega <- omega + (1 - h / (lag + 1)) * (g + t(g))
  }
  omega
}

sandwich_se <- function(score, hessian = NULL, h_inv = NULL,
                        hac_lag = 0, coef_names = NULL) {
  score <- as.matrix(score)
  if (is.null(h_inv)) {
    if (is.null(hessian)) stop("provide hessian or h_inv")
    h_inv <- safe_inv(hessian)
  }
  omega <- if (hac_lag > 0) meat_hac(score, hac_lag) else meat_iid(score)
  vcov_sqrt_t <- sym(h_inv %*% omega %*% h_inv)
  vcov <- vcov_sqrt_t / nrow(score)
  se <- sqrt(pmax(diag(vcov), 0))
  if (!is.null(coef_names)) names(se) <- coef_names
  list(se = se, vcov = vcov, vcov_sqrt_t = vcov_sqrt_t,
       omega = omega, h_inv = h_inv)
}

expectile_sandwich <- function(y, x, beta, tau, hac_lag = 0, coef_names = NULL) {
  x <- as.matrix(x)
  resid <- as.numeric(y - x %*% beta)
  score <- x * as.numeric(rho_dot(resid, tau))
  hessian <- crossprod(x, x * omega_tau(resid, tau)) / nrow(x)
  sandwich_se(score, hessian = hessian, hac_lag = hac_lag, coef_names = coef_names)
}

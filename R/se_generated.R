as_x_list <- function(x, k) {
  if (is.list(x)) return(x)
  replicate(k, x, simplify = FALSE)
}

stage1_influence <- function(x_list, e_hat, tau) {
  # First-stage influence used by the generated-residual correction:
  # L_{j,t}^{(1)} = x_{j,t} H_j^{-1} w_tau(e_{j,t}) e_{j,t},
  # H_j = T^{-1} sum_t w_tau(e_{j,t}) x_{j,t} x_{j,t}'.
  k <- ncol(e_hat)
  lapply(seq_len(k), function(j) {
    xj <- x_list[[j]]
    wj <- expectile_weight(e_hat[, j], tau)
    wj_inv <- safe_inv(crossprod(xj, xj * wj) / length(wj))
    (xj %*% wj_inv) * as.numeric(wj * e_hat[, j])
  })
}

generated_se_i <- function(i, x, e_hat, eta_i, theta_i, tau,
                           hac_lag = 0, coef_names = NULL, residual_sources = NULL,
                           first_stage_x = NULL) {
  # Generated-residual covariance for equation i.
  # Moment: 0 = T^{-1} sum_t z_{i,t} psi_tau(eta_{i,t}),
  # z_{i,t} = (x_{i,t}', ehat_{j,t}: j in J_i)'.
  #
  # Paper expansion implemented below:
  # sqrt(T)(theta_hat_i - theta_i)
  #   = H_i^{-1} T^{-1/2} sum_t L_{i,t} + o_p(1),
  # H_i = T^{-1} sum_t omega_tau(eta_{i,t}) z_{i,t} z_{i,t}'.
  #
  # L_{i,t} = z_{i,t} psi_tau(eta_{i,t})
  #           + sum_{j in J_i} L_{j,t}^{(1)} B_{ij}',
  # B_{ij} = delta_{ij} T^{-1} sum_t omega_tau(eta_{i,t}) z_{i,t} x_{j,t}'
  #          - 1{coordinate = resid_j} T^{-1} sum_t psi_tau(eta_{i,t}) x_{j,t}'.
  #
  # Var(theta_hat_i) = H_i^{-1} Omega_i H_i^{-1} / T,
  # Omega_i = Var(T^{-1/2} sum_t L_{i,t}), estimated by IID or HAC meat.
  e_hat <- as.matrix(e_hat)
  k <- ncol(e_hat)
  x_list <- as_x_list(x, k)
  first_stage_x_list <- as_x_list(first_stage_x %||% x, k)
  other <- residual_sources %||% setdiff(seq_len(k), i)
  other <- setdiff(other, i)
  xi <- x_list[[i]]
  p <- ncol(xi)
  n <- nrow(xi)

  z <- cbind(xi, e_hat[, other, drop = FALSE])
  psi <- rho_dot(eta_i, tau)
  h_inv <- safe_inv(crossprod(z, z * omega_tau(eta_i, tau)) / n)
  base_l <- z * as.numeric(psi)
  l1 <- stage1_influence(first_stage_x_list, e_hat, tau)
  gen_l <- matrix(0, n, ncol(z))

  for (r in seq_along(other)) {
    j <- other[r]
    xj <- first_stage_x_list[[j]]
    delta_j <- theta_i[p + r]
    b_j <- delta_j * crossprod(z, xj * omega_tau(eta_i, tau)) / n
    b_j[p + r, ] <- b_j[p + r, ] - colMeans(xj * as.numeric(psi))
    gen_l <- gen_l + l1[[j]] %*% t(b_j)
  }

  l_fg <- base_l + gen_l
  omega <- if (hac_lag > 0) meat_hac(l_fg, hac_lag) else meat_iid(l_fg)
  vcov_sqrt_t <- sym(h_inv %*% omega %*% h_inv)
  vcov <- vcov_sqrt_t / n
  se <- sqrt(pmax(diag(vcov), 0))
  if (!is.null(coef_names)) names(se) <- coef_names
  list(se = se, vcov = vcov, vcov_sqrt_t = vcov_sqrt_t,
       influence = l_fg, base_influence = base_l,
       generated_influence = gen_l, h_inv = h_inv, omega = omega)
}

generated_se <- function(x, e_hat, eta_hat, theta_hat, tau,
                         hac_lag = 0, coef_names = NULL, residual_sources = NULL,
                         first_stage_x = NULL) {
  e_hat <- as.matrix(e_hat)
  eta_hat <- as.matrix(eta_hat)
  k <- ncol(e_hat)
  x_list <- as_x_list(x, k)
  theta_list <- if (is.list(theta_hat)) {
    theta_hat
  } else if (is.matrix(theta_hat) && nrow(theta_hat) == k) {
    lapply(seq_len(k), function(i) as.numeric(theta_hat[i, ]))
  } else if (is.matrix(theta_hat) && ncol(theta_hat) == k) {
    lapply(seq_len(k), function(i) as.numeric(theta_hat[, i]))
  } else {
    stop("theta_hat must be a list or a matrix with one equation per row/column")
  }
  lapply(seq_len(k), function(i) {
    names_i <- if (is.list(coef_names)) coef_names[[i]] else coef_names
    src_i <- if (is.list(residual_sources)) residual_sources[[i]] else residual_sources
    generated_se_i(i, x_list, e_hat, eta_hat[, i], theta_list[[i]],
                   tau, hac_lag = hac_lag, coef_names = names_i,
                   residual_sources = src_i, first_stage_x = first_stage_x)
  })
}

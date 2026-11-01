sals <- function(y, spec = "ABSQ", tau = 0.99855, q = 1,
                 own_only = FALSE, hac_lag = 0) {
  design <- see_build_design(y, spec = spec, q = q)
  yy <- design$y
  equations <- colnames(yy)
  k <- ncol(yy)

  coef_rows <- gof_rows <- diag_rows <- vector("list", k)
  fits <- x_list <- vcov_list <- term_names <- term_maps <- vector("list", k)
  fitted_mat <- resid_mat <- matrix(NA_real_, nrow(yy), k, dimnames = dimnames(yy))

  for (i in seq_len(k)) {
    keep <- if (own_only) is.na(design$colmap$source_index) | design$colmap$source_index == i else rep(TRUE, ncol(design$x))
    x_i <- design$x[, keep, drop = FALSE]
    map_i <- design$colmap[keep, , drop = FALSE]
    y_i <- yy[, i]
    check_fixed_design(x_i, paste0("SALS equation ", equations[i]))
    fit <- fit_expectile(y_i, x_i, tau)
    sw <- expectile_sandwich(y_i, x_i, fit$beta, tau, hac_lag, colnames(x_i))

    fitted_mat[, i] <- fit$fitted
    resid_mat[, i] <- fit$resid
    fits[[i]] <- fit
    x_list[[i]] <- x_i
    vcov_list[[i]] <- sw$vcov
    term_names[[i]] <- colnames(x_i)
    term_maps[[i]] <- map_i

    coef_rows[[i]] <- data.frame(
      stage = "SALS", equation = equations[i], term = colnames(x_i),
      term_type = ifelse(map_i$transform == "intercept", "intercept", "lag"),
      source_equation = map_i$source_equation,
      estimate = fit$beta, sandwich_se = sw$se,
      t = fit$beta / sw$se,
      p = 2 * stats::pnorm(abs(fit$beta / sw$se), lower.tail = FALSE),
      iter = fit$iter, converged = fit$converged, score_max = fit$score_max,
      row.names = NULL
    )

    gof_rows[[i]] <- cbind(
      stage = "SALS", equation = equations[i], tau = tau, spec = design$spec,
      see_gof(y_i, fit$fitted, ncol(x_i), tau,
              p_effective = qr(x_i)$rank)
    )
    diag_rows[[i]] <- cbind(stage = "SALS", equation = equations[i], residual_diagnostics(fit$resid))
  }

  structure(
    list(
      call = match.call(), tau = tau, spec = design$spec, q = q,
      own_only = own_only, y = yy, dates = design$dates, design = design,
      x_list = x_list, fits = fits, fitted = fitted_mat, resid = resid_mat,
      coef = do.call(rbind, coef_rows), gof = do.call(rbind, gof_rows),
      diagnostics = do.call(rbind, diag_rows),
      vcov_sandwich = vcov_list, term_names = term_names, term_maps = term_maps
    ),
    class = "expectsee_sals"
  )
}

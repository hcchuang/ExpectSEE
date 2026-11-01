sfgals <- function(stage1, se = c("generated", "sandwich", "none"), hac_lag = 0) {
  se <- match.arg(se)
  yy <- stage1$y
  e_hat <- stage1$resid
  equations <- colnames(yy)
  k <- ncol(yy)

  coef_rows <- gof_rows <- diag_rows <- vector("list", k)
  fits <- x_list <- vcov_sand <- vcov_gen <- term_names <- vector("list", k)
  fitted_mat <- resid_mat <- matrix(NA_real_, nrow(yy), k, dimnames = dimnames(yy))

  for (i in seq_len(k)) {
    other <- setdiff(seq_len(k), i)
    z_i <- cbind(stage1$x_list[[i]], e_hat[, other, drop = FALSE])
    res_terms <- paste0("resid_", equations[other])
    colnames(z_i) <- c(stage1$term_names[[i]], res_terms)
    y_i <- yy[, i]
    check_fixed_design(z_i, paste0("SFGALS equation ", equations[i]))
    fit <- fit_expectile(y_i, z_i, stage1$tau)
    sw <- expectile_sandwich(y_i, z_i, fit$beta, stage1$tau, hac_lag, colnames(z_i))
    gen <- generated_se_i(i, stage1$x_list, e_hat, fit$resid, fit$beta,
                          stage1$tau, hac_lag, colnames(z_i))

    fitted_mat[, i] <- fit$fitted
    resid_mat[, i] <- fit$resid
    fits[[i]] <- fit
    x_list[[i]] <- z_i
    vcov_sand[[i]] <- sw$vcov
    vcov_gen[[i]] <- gen$vcov
    term_names[[i]] <- colnames(z_i)

    se_use <- if (se == "sandwich") sw$se else gen$se
    coef_rows[[i]] <- data.frame(
      stage = "SFGALS", equation = equations[i], term = colnames(z_i),
      term_type = c(
        ifelse(stage1$term_maps[[i]]$transform == "intercept", "intercept", "lag"),
        rep("generated_residual", length(other))
      ),
      source_equation = c(stage1$term_maps[[i]]$source_equation, equations[other]),
      estimate = fit$beta,
      sandwich_se = sw$se,
      generated_se = gen$se,
      t = fit$beta / se_use,
      p = 2 * stats::pnorm(abs(fit$beta / se_use), lower.tail = FALSE),
      iter = fit$iter, converged = fit$converged, score_max = fit$score_max,
      row.names = NULL
    )

    gof_rows[[i]] <- cbind(
      stage = "SFGALS", equation = equations[i], tau = stage1$tau, spec = stage1$spec,
      see_gof(y_i, fit$fitted, ncol(z_i), stage1$tau,
              p_effective = qr(z_i)$rank)
    )
    diag_rows[[i]] <- cbind(stage = "SFGALS", equation = equations[i], residual_diagnostics(fit$resid))
  }

  structure(
    list(
      call = match.call(), tau = stage1$tau, spec = stage1$spec, q = stage1$q,
      first_stage = stage1, y = yy, dates = stage1$dates, fits = fits,
      x_list = x_list,
      fitted = fitted_mat, resid = resid_mat, coef = do.call(rbind, coef_rows),
      gof = do.call(rbind, gof_rows), diagnostics = do.call(rbind, diag_rows),
      vcov_sandwich = vcov_sand, vcov_generated = vcov_gen,
      term_names = term_names, se = se
    ),
    class = "expectsee_sfgals"
  )
}

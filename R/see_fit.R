see_fit <- function(y, spec = "ABSQ", tau = 0.99855, q = 1,
                    estimator = c("sfgals", "sals", "hdsfgals", "hdsals"),
                    se = c("generated", "sandwich", "none"),
                    own_only = FALSE, hac_lag = 0,
                    coefficient_fit = c("lasso", "post_lasso"),
                    nested_residual_cv = NULL,
                    stage2_lag_candidates = c("full", "stage1_selected"),
                    stage2_penalty = c("all", "stage1_lags_unpenalized"), ...) {
  estimator <- match.arg(estimator)
  se <- match.arg(se)
  coefficient_fit <- match.arg(coefficient_fit)
  stage2_lag_candidates <- match.arg(stage2_lag_candidates)
  stage2_penalty <- match.arg(stage2_penalty)
  if (estimator %in% c("hdsals", "hdsfgals")) {
    stage1 <- hdsals(
      y, spec = spec, tau = tau, q = q, own_only = own_only,
      coefficient_fit = coefficient_fit, hac_lag = hac_lag, ...
    )
    stage2 <- if (estimator == "hdsfgals") hdsfgals(
      stage1, coefficient_fit = coefficient_fit,
      se = if (coefficient_fit == "lasso") "none" else se,
      hac_lag = hac_lag, nested_residual_cv = nested_residual_cv,
      stage2_lag_candidates = stage2_lag_candidates,
      stage2_penalty = stage2_penalty, ...
    ) else NULL
    active <- stage2 %||% stage1
  } else {
    stage1 <- sals(y, spec = spec, tau = tau, q = q, own_only = own_only, hac_lag = hac_lag)
    stage2 <- if (estimator == "sfgals") sfgals(stage1, se = se, hac_lag = hac_lag) else NULL
    active <- stage2 %||% stage1
  }

  structure(
    list(
      call = match.call(), tau = tau, spec = stage1$spec, q = q,
      estimator = estimator,
      coefficient_fit = if (estimator %in% c("hdsals", "hdsfgals")) coefficient_fit else NA_character_,
      coefficient_fit_label = if (estimator %in% c("hdsals", "hdsfgals")) active$coefficient_fit_label else NA_character_,
      se = if (estimator %in% c("hdsals", "hdsfgals")) active$se else se,
      sals = if (estimator %in% c("sals", "sfgals")) stage1 else NULL,
      sfgals = if (estimator == "sfgals") stage2 else NULL,
      hdsals = if (estimator %in% c("hdsals", "hdsfgals")) stage1 else NULL,
      hdsfgals = if (estimator == "hdsfgals") stage2 else NULL,
      coef = active$coef, gof = active$gof, diagnostics = active$diagnostics,
      fitted = active$fitted, resid = active$resid,
      support = active$support, edges = active$edges,
      lambda_path = active$lambda_path
    ),
    class = "expectsee"
  )
}

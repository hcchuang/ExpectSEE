see_aicr_simple <- function(loss_sum, n, p_effective, reference_scale) {
  if (any(!is.finite(c(loss_sum, n, p_effective, reference_scale))) ||
      any(n <= 0) || any(p_effective < 0) || any(reference_scale <= 0)) {
    stop("AICR inputs must be finite, with n and reference_scale positive")
  }
  2 * loss_sum / (n * reference_scale) + 2 * p_effective / n
}

see_aicr_reference_scale <- function(resid, tau) {
  resid <- as.matrix(resid)
  vapply(seq_len(ncol(resid)), function(i) {
    e <- resid[, i]
    ok <- is.finite(e)
    if (!any(ok)) return(NA_real_)
    sum(expectile_loss(e[ok], tau)) / sum(expectile_weight(e[ok], tau))
  }, numeric(1))
}

see_gof <- function(y, yhat, p, tau, reference_scale = NULL, p_effective = p) {
  y <- as.numeric(y)
  yhat <- as.numeric(yhat)
  if (!length(y) || length(y) != length(yhat)) {
    stop("y and yhat must be nonempty and have the same length")
  }
  if (length(p) != 1 || !is.finite(p) || p < 0) {
    stop("p must be one nonnegative finite number, including the intercept")
  }
  if (length(p_effective) != 1 || !is.finite(p_effective) || p_effective < 0) {
    stop("p_effective must be one nonnegative finite number")
  }
  if (!is.null(reference_scale) && length(reference_scale) != 1) {
    stop("reference_scale must be one number for one equation")
  }
  resid <- y - yhat
  loss <- expectile_loss(resid, tau)
  n <- length(y)
  null_loss <- sum(fit_expectile(y, matrix(1, n, 1), tau)$loss)
  null_tol <- .Machine$double.eps * max(1, sum(y^2))
  r2 <- if (null_loss > null_tol) 1 - sum(loss) / null_loss else NA_real_
  adj_r2 <- if (n > p && is.finite(r2)) 1 - (1 - r2) * (n - 1) / (n - p) else NA_real_
  aicr <- if (is.null(reference_scale)) NA_real_ else {
    see_aicr_simple(sum(loss), n, p_effective, reference_scale)
  }
  data.frame(
    n = n, p = p, p_effective = p_effective,
    loss_sum = sum(loss),
    loss_mean = mean(loss),
    R2 = r2,
    adjR2 = adj_r2,
    reference_scale = reference_scale %||% NA_real_,
    AICR_simple = aicr,
    AICR_version = "simple_v1",
    stringsAsFactors = FALSE
  )
}

residual_diagnostics <- function(resid) {
  resid <- as.numeric(resid)
  data.frame(
    mean = mean(resid),
    sd = stats::sd(resid),
    min = min(resid),
    max = max(resid),
    p05 = as.numeric(stats::quantile(resid, 0.05)),
    p50 = as.numeric(stats::quantile(resid, 0.50)),
    p95 = as.numeric(stats::quantile(resid, 0.95)),
    acf1 = if (length(resid) > 2) stats::acf(resid, plot = FALSE, lag.max = 1)$acf[2] else NA_real_,
    stringsAsFactors = FALSE
  )
}

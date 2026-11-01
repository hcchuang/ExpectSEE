hac_lrv <- function(x, lag = NULL) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  z <- x - mean(x)
  if (is.null(lag)) lag <- floor(4 * (length(z) / 100)^(2 / 9))
  out <- sum(z * z) / length(z)
  for (h in seq_len(lag)) {
    out <- out + 2 * (1 - h / (lag + 1)) * sum(z[(h + 1):length(z)] * z[1:(length(z) - h)]) / length(z)
  }
  max(out, 0)
}

see_parzen_lrv <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  z <- x - mean(x)
  m <- min(ceiling(2 * sqrt(length(z))), length(z) - 1)
  gam <- as.numeric(stats::acf(z, lag.max = m, type = "covariance",
                               plot = FALSE, na.action = stats::na.fail)$acf)
  u <- seq_len(m) / m
  kernel <- ifelse(u <= 0.5, 1 - 6 * u^2 + 6 * u^3, 2 * (1 - u)^3)
  lam <- c(1, 2 * kernel)
  max(as.numeric(gam %*% lam), 0)
}

see_nz_score <- function(y, yhat, tau, h = 1, positive_guard = FALSE) {
  check_tau(tau)
  h <- match.arg(as.character(h), c("1", "0"))
  x <- as.numeric(y)
  r <- as.numeric(yhat)
  if (length(x) != length(r)) stop("y and yhat must have the same length")
  s <- rep(NA_real_, length(r))
  ok <- is.finite(x) & is.finite(r)

  if (h == "1" && !positive_guard) {
    s[ok] <- (1 - tau) * r[ok] * (r[ok] - 2 * x[ok]) -
      as.numeric(x[ok] > r[ok]) * (1 - 2 * tau) * (x[ok] - r[ok])^2
    return(s)
  }

  s[ok] <- 0
  ok <- ok & r > 0
  ind1 <- ok & x > r
  ind2 <- ok & x < r
  if (h == "1") {
    s[ind1] <- (-(1 - 2 * tau) * (x - r)^2 + (1 - tau) * r * (r - 2 * x))[ind1]
    s[ind2] <- ((1 - tau) * r * (r - 2 * x))[ind2]
  } else {
    s[ind1] <- ((1 - 2 * tau) * (log(x / r) + 1 - x / r) +
      (1 - tau) * (log(r) - 1 + x / r))[ind1]
    s[ind2] <- ((1 - tau) * (log(r) - 1 + x / r))[ind2]
  }
  s
}

score_test <- function(diff, lrv, alpha, mean_name) {
  diff <- as.numeric(diff)
  diff <- diff[is.finite(diff)]
  mean_diff <- mean(diff)
  if (is.finite(lrv) && lrv > 0) {
    z <- mean_diff / sqrt(lrv / length(diff))
  } else if (is.finite(lrv) && lrv == 0 && abs(mean_diff) < 1e-14) {
    z <- 0
  } else if (is.finite(lrv) && lrv == 0) {
    z <- sign(mean_diff) * Inf
  } else {
    z <- NA_real_
  }
  p_better <- stats::pnorm(z)
  p_worse <- stats::pnorm(z, lower.tail = FALSE)
  out <- data.frame(
    n_obs = length(diff), z = z,
    p_internal_better = p_better,
    p_internal_worse = p_worse,
    decision = traffic_light(p_better, p_worse, alpha),
    row.names = NULL
  )
  out[[mean_name]] <- mean_diff
  out[c("n_obs", mean_name, "z", "p_internal_better", "p_internal_worse", "decision")]
}

dm_test_loss <- function(internal_loss, standard_loss, lag = NULL, alpha = 0.05) {
  diff <- as.numeric(internal_loss) - as.numeric(standard_loss)
  diff <- diff[is.finite(diff)]
  score_test(diff, hac_lrv(diff, lag), alpha, "mean_loss_diff")
}

see_nz_tlm_test <- function(internal_yhat, standard_yhat, y, tau, alpha = 0.05,
                            h = 1, positive_guard = FALSE) {
  s1 <- see_nz_score(y, internal_yhat, tau, h, positive_guard)
  s0 <- see_nz_score(y, standard_yhat, tau, h, positive_guard)
  out <- score_test(s1 - s0, see_parzen_lrv(s1 - s0), alpha, "mean_score_diff")
  out$p_H0_plus <- out$p_internal_better
  out$p_H0_minus <- out$p_internal_worse
  out
}

see_traffic_light <- function(p_better, p_worse, alpha = 0.05) {
  ifelse(p_better <= alpha, "green", ifelse(p_worse <= alpha, "red", "yellow"))
}

traffic_light <- see_traffic_light

see_compare_forecasts <- function(forecast, tau, alpha = 0.05, lag = NULL,
                                  equation = "equation", model = "model",
                                  y = "y", yhat = "yhat", date = "date",
                                  method = c("nz", "dm"), h = 1,
                                  positive_guard = FALSE) {
  method <- match.arg(method)
  need <- c(equation, model, y, yhat)
  if (!all(need %in% names(forecast))) stop("forecast is missing required columns")
  forecast$loss <- expectile_loss(forecast[[y]] - forecast[[yhat]], tau)
  if (!date %in% names(forecast)) forecast[[date]] <- seq_len(nrow(forecast))

  rows <- list()
  by_eq <- split(forecast, forecast[[equation]])
  for (eq in names(by_eq)) {
    g <- by_eq[[eq]]
    models <- unique(g[[model]])
    for (m1 in models) for (m0 in models) {
      if (identical(m1, m0)) next
      a <- g[g[[model]] == m1, c(date, y, yhat, "loss")]
      b <- g[g[[model]] == m0, c(date, yhat, "loss")]
      names(a) <- c("date", "y", "yhat_internal", "loss_internal")
      names(b) <- c("date", "yhat_standard", "loss_standard")
      common <- merge(a, b, by = "date", suffixes = c("_internal", "_standard"))
      stat <- if (method == "nz") {
        see_nz_tlm_test(common$yhat_internal, common$yhat_standard, common$y,
                        tau, alpha, h, positive_guard)
      } else {
        dm_test_loss(common$loss_internal, common$loss_standard, lag, alpha)
      }
      rows[[length(rows) + 1]] <- cbind(
        equation = eq, traffic_method = method,
        internal_model = m1, standard_model = m0, stat
      )
    }
  }
  do.call(rbind, rows)
}

see_backtest <- see_compare_forecasts

see_reference_scale <- function(fit, tau) {
  resid <- as.matrix(fit$resid)
  equations <- colnames(resid) %||% as.character(fit$gof$equation)
  reference_scale <- see_aicr_reference_scale(resid, tau)
  data.frame(equation = equations, reference_scale = reference_scale,
             stringsAsFactors = FALSE)
}

score_aicr <- function(gof, reference) {
  reference_scale <- reference$reference_scale[match(gof$equation, reference$equation)]
  if (any(!is.finite(reference_scale) | reference_scale <= 0)) {
    stop("ABSQ reference scale is missing or nonpositive for at least one equation")
  }
  gof$reference_scale <- reference_scale
  gof$AICR_simple <- see_aicr_simple(
    gof$loss_sum, gof$n, gof$p_effective, reference_scale
  )
  gof$AICR_version <- "simple_v1"
  gof
}

see_model_selection <- function(y, specs = see_specs(), taus = 0.99855,
                                estimator = c("sfgals", "sals", "hdsfgals", "hdsals"),
                                criterion = c("AICR_simple", "loss_mean", "R2", "adjR2"),
                                q = 1, own_only = FALSE, se = "generated", ...) {
  estimator <- match.arg(estimator)
  criterion <- match.arg(criterion)
  specs <- unique(vapply(specs, normalize_spec, character(1)))
  grid <- do.call(rbind, lapply(specs, function(sp) {
    data.frame(spec = sp, tau = taus, stringsAsFactors = FALSE)
  }))

  # Fit first. AICR is comparison-dependent and is scored only after the
  # common ABSQ reference scale is available for each tau.
  fits <- lapply(seq_len(nrow(grid)), function(i) {
    see_fit(y, spec = grid$spec[i], tau = grid$tau[i], q = q,
            estimator = estimator, own_only = own_only, se = se, ...)
  })
  references <- vector("list", length(unique(taus)))
  names(references) <- as.character(unique(taus))

  for (ta in unique(taus)) {
    id <- which(grid$tau == ta)
    ref_id <- id[grid$spec[id] == "ABSQ"]
    ref_fit <- if (length(ref_id)) fits[[ref_id[1]]] else {
      see_fit(y, spec = "ABSQ", tau = ta, q = q, estimator = estimator,
              own_only = own_only, se = se, ...)
    }
    reference <- see_reference_scale(ref_fit, ta)
    reference$tau <- ta
    reference$reference_spec <- "ABSQ"
    reference$reference_estimator <- estimator
    reference$reference_coefficient_fit <- ref_fit$coefficient_fit %||% NA_character_
    references[[as.character(ta)]] <- reference
    for (i in id) {
      fits[[i]]$gof <- score_aicr(fits[[i]]$gof, reference)
    }
  }

  finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    g <- fits[[i]]$gof
    data.frame(
      spec = grid$spec[i], tau = grid$tau[i], estimator = estimator,
      coefficient_fit = fits[[i]]$coefficient_fit %||% NA_character_,
      criterion = criterion, score = finite_mean(g[[criterion]]),
      loss_mean = finite_mean(g$loss_mean), R2 = finite_mean(g$R2),
      adjR2 = finite_mean(g$adjR2), AICR_simple = finite_mean(g$AICR_simple),
      AICR_version = "simple_v1",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  equation_gof <- do.call(rbind, lapply(fits, function(fit) fit$gof))
  valid <- which(is.finite(out$score))
  if (!length(valid)) stop("No finite model-selection scores were produced")
  pick <- if (criterion %in% c("R2", "adjR2")) which.max else which.min
  best <- do.call(rbind, lapply(unique(out$tau), function(ta) {
    id <- which(out$tau == ta & is.finite(out$score))
    if (!length(id)) return(NULL)
    out[id[pick(out$score[id])], , drop = FALSE]
  }))
  structure(
    list(results = out, best = best,
         equation_gof = equation_gof,
         reference_scale = do.call(rbind, references)),
    class = "expectsee_selection"
  )
}

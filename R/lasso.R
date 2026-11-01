need_sales <- function() {
  if (!requireNamespace("SALES", quietly = TRUE)) {
    stop("HDSALS/HDSFGALS require SALES for expectile lasso.", call. = FALSE)
  }
  invisible(TRUE)
}

sales_call <- function(name) {
  need_sales()
  utils::getFromNamespace(name, "SALES")
}

see_lambda_grid <- function(nlambda = 100, lambda_max = 1, lambda_min = 1e-10) {
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (!is.finite(lambda_max) || !is.finite(lambda_min) ||
      lambda_max <= lambda_min || lambda_min <= 0) {
    stop("lambda grid requires 0 < lambda_min < lambda_max")
  }
  exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
}

see_ts_splits <- function(n, initial = floor(n / 2), assess = max(1, floor(n / 10)),
                          step = assess, type = c("fixed_window", "expanding_window")) {
  type <- match.arg(type)
  if (initial < 3 || assess < 1 || initial + assess > n) stop("invalid time-series CV window sizes")
  stops <- seq(initial, n - assess, by = step)
  lapply(stops, function(last_train) {
    train <- if (type == "expanding_window") seq_len(last_train) else (last_train - initial + 1):last_train
    list(train = train, test = (last_train + 1):(last_train + assess))
  })
}

safe_ts_splits <- function(n, type, initial = NULL, assess = NULL, step = NULL) {
  assess <- assess %||% max(1, floor(n / 10))
  starts <- unique(c(initial %||% floor(n / 2), floor(0.6 * n), n - 2 * assess))
  for (init in starts) {
    init <- max(3, min(as.integer(init), n - assess - 1))
    ass <- max(1, min(as.integer(assess), n - init))
    st <- step %||% ass
    out <- tryCatch(see_ts_splits(n, init, ass, st, type), error = function(e) NULL)
    if (length(out)) return(out)
  }
  stop("could not build nested time-series CV splits")
}

choose_lambda <- function(lambda, cvm, cvsd) {
  lambda_min <- max(lambda[cvm <= min(cvm, na.rm = TRUE)])
  sd_min <- cvsd[match(lambda_min, lambda)]
  if (!is.finite(sd_min)) sd_min <- 0  # one validation split: no standard error, one-SE = minimum
  se_cut <- cvm[match(lambda_min, lambda)] + sd_min
  lambda_1se <- max(lambda[cvm <= se_cut])
  list(lambda.min = lambda_min, lambda.1se = lambda_1se)
}

lasso_nzero <- function(fit, lambda) {
  nz <- tryCatch(lengths(stats::coef(fit, type = "nonzero")), error = function(e) integer())
  if (length(nz) == length(lambda)) return(as.integer(nz))
  cf <- tryCatch(as.matrix(stats::coef(fit, s = lambda)), error = function(e) NULL)
  if (is.null(cf)) return(rep(NA_integer_, length(lambda)))
  keep <- rownames(cf) != "(Intercept)"
  as.integer(colSums(abs(cf[keep, , drop = FALSE]) > 0))
}

cv_expectile_lasso_ts <- function(y, x, tau, splits, penalty_factor, lambda = NULL,
                                  lambda2 = 0, standardize = TRUE, intercept = TRUE, ...) {
  if (length(splits) < 1) stop("time-series CV needs at least one split")
  ernet <- sales_call("ernet")

  full_fit <- ernet(x = x, y = y, tau = tau, pf = penalty_factor, lambda = lambda,
                    lambda2 = lambda2, standardize = standardize, intercept = intercept, ...)
  lambda <- full_fit$lambda
  split_loss <- matrix(NA_real_, length(splits), length(lambda))

  for (s in seq_along(splits)) {
    id <- splits[[s]]
    # Forward validation: train only on past observations and score the next block.
    fit_s <- ernet(x = x[id$train, , drop = FALSE], y = y[id$train], tau = tau,
                   pf = penalty_factor, lambda = lambda, lambda2 = lambda2,
                   standardize = standardize, intercept = intercept, ...)
    resid <- y[id$test] - stats::predict(fit_s, newx = x[id$test, , drop = FALSE], s = lambda)
    split_loss[s, ] <- colMeans(abs(tau - (resid < 0)) * resid^2)
  }

  cvm <- colMeans(split_loss, na.rm = TRUE)
  cvsd <- apply(split_loss, 2, stats::sd, na.rm = TRUE) / sqrt(nrow(split_loss))
  picked <- choose_lambda(lambda, cvm, cvsd)
  nzero <- lasso_nzero(full_fit, lambda)
  structure(
    c(
      list(lambda = lambda, cvm = cvm, cvsd = cvsd, nzero = nzero,
           ernet.fit = full_fit, splits = splits, cv_method = "time_series"),
      picked
    ),
    class = "expectsee_ts_cv_ernet"
  )
}

predict_lasso_beta <- function(beta, x) {
  b0 <- if ("(Intercept)" %in% names(beta)) beta["(Intercept)"] else if ("const" %in% names(beta)) beta["const"] else 0
  terms <- intersect(setdiff(names(beta), c("(Intercept)", "const")), colnames(x))
  if (!length(terms)) return(rep(as.numeric(b0), nrow(x)))
  as.numeric(b0 + x[, terms, drop = FALSE] %*% beta[terms])
}

nested_stage1_for_split <- function(stage1, fit_idx, pred_idx, cv, lambda_choice,
                                    nfolds, standardize, lambda, lambda2,
                                    zero_tol, initial, assess, step,
                                    coefficient_fit = "lasso", ...) {
  yy <- stage1$y
  equations <- colnames(yy)
  out <- matrix(NA_real_, nrow(yy), ncol(yy), dimnames = dimnames(yy))
  selected <- stats::setNames(vector("list", ncol(yy)), equations)
  rows <- path_rows <- vector("list", ncol(yy))

  for (j in seq_along(equations)) {
    keep <- if (stage1$own_only) stage1$design$colmap$source_index == j else rep(TRUE, ncol(stage1$x_raw))
    x_j <- stage1$x_raw[, keep, drop = FALSE]
    y_j <- yy[, j]
    inner <- safe_ts_splits(length(fit_idx), cv, initial, assess, step)
    lambda_j <- if (is.list(lambda)) lambda[[j]] else lambda
    lasso <- fit_expectile_lasso(
      y_j[fit_idx], x_j[fit_idx, , drop = FALSE], stage1$tau, lambda_choice,
      nfolds = nfolds, standardize = standardize, lambda = lambda_j,
      lambda2 = lambda2, zero_tol = zero_tol, cv = cv, ts_splits = inner, ...
    )
    selected[[j]] <- intersect(lasso$selected_terms, colnames(x_j))
    coefficient <- selected_coefficient_fit(
      lasso, y_j[fit_idx], x_j[fit_idx, , drop = FALSE], selected[[j]],
      stage1$tau, coefficient_fit
    )
    pred <- predict_lasso_beta(coefficient$fit$beta, x_j)
    use <- c(fit_idx, pred_idx)
    out[use, j] <- y_j[use] - pred[use]
    rows[[j]] <- data.frame(
      stage1_equation = equations[j], lambda = lasso$lambda,
      n_selected = length(selected[[j]]), n_train = length(fit_idx),
      coefficient_fit = coefficient_fit,
      pred_start = min(pred_idx), pred_end = max(pred_idx),
      row.names = NULL
    )
    path_rows[[j]] <- cbind(
      stage1_equation = equations[j], n_train = length(fit_idx),
      pred_start = min(pred_idx), pred_end = max(pred_idx),
      lasso$path
    )
  }
  list(
    resid = out, selected = selected, info = do.call(rbind, rows),
    path = do.call(rbind, path_rows)
  )
}

stage2_lag_penalty <- function(lag_terms, forced_lags = character()) {
  if (!length(lag_terms)) return(numeric())
  as.numeric(!lag_terms %in% forced_lags)
}

stage2_design <- function(stage1, i, residuals, lag_terms = NULL,
                          forced_lags = character(), residual_penalty = 1) {
  equations <- colnames(stage1$y)
  other <- setdiff(seq_along(equations), i)
  res_other <- residuals[, other, drop = FALSE]
  colnames(res_other) <- paste0("resid_", equations[other])
  lag_terms <- lag_terms %||% colnames(stage1$x_raw)
  missing <- setdiff(lag_terms, colnames(stage1$x_raw))
  if (length(missing)) {
    stop("Stage 2 lag candidates not found in Stage 1 design: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  x_lag <- stage1$x_raw[, lag_terms, drop = FALSE]
  lag_pf <- stage2_lag_penalty(lag_terms, forced_lags)
  list(
    x = cbind(x_lag, res_other),
    pf = c(lag_pf, rep(residual_penalty, ncol(res_other))),
    lag_terms = lag_terms,
    forced_lags = intersect(lag_terms, forced_lags),
    residual_terms = colnames(res_other),
    residual_sources = other
  )
}

fit_hdsfgals_lasso_nested <- function(stage1, i, lambda_choice, standardize,
                                      lambda, lambda2, zero_tol, cv, ts_splits,
                                      initial, assess, step, nfolds,
                                      stage2_lag_candidates, lag_terms,
                                      forced_lags, residual_penalty,
                                      coefficient_fit = "lasso", ...) {
  if (cv == "iid_kfold") stop("nested_residual_cv requires fixed_window or expanding_window CV")
  ernet <- sales_call("ernet")
  yy <- stage1$y
  y_i <- yy[, i]
  splits <- ts_splits %||% safe_ts_splits(length(y_i), cv, initial, assess, step)
  full_d <- stage2_design(stage1, i, stage1$resid, lag_terms, forced_lags, residual_penalty)
  full_fit <- ernet(x = full_d$x, y = y_i, tau = stage1$tau, pf = full_d$pf,
                    lambda = lambda, lambda2 = lambda2,
                    standardize = standardize, intercept = TRUE, ...)
  lambda <- full_fit$lambda
  stage1_tuning <- stage1$tuning %||% list()
  stage1_lambda <- lapply(stage1$selection_fits %||% stage1$fits,
                          function(z) z$path$lambda)
  loss <- matrix(NA_real_, length(splits), length(lambda))
  fold_rows <- nested_rows <- nested_path_rows <- vector("list", length(splits))

  for (s in seq_along(splits)) {
    id <- splits[[s]]
    nested <- nested_stage1_for_split(
      stage1, id$train, id$test, cv,
      stage1_tuning$lambda_choice %||% lambda_choice,
      stage1_tuning$nfolds %||% nfolds,
      stage1_tuning$standardize %||% standardize,
      stage1_lambda, stage1_tuning$lambda2 %||% 0,
      stage1_tuning$zero_tol %||% zero_tol,
      stage1_tuning$initial, stage1_tuning$assess,
      stage1_tuning$step, coefficient_fit, ...
    )
    fold_lag_terms <- if (stage2_lag_candidates == "stage1_selected") nested$selected[[i]] else lag_terms
    fold_forced_lags <- if (length(forced_lags)) nested$selected[[i]] else character()
    d <- stage2_design(stage1, i, nested$resid, fold_lag_terms,
                       fold_forced_lags, residual_penalty)
    fit_s <- ernet(x = d$x[id$train, , drop = FALSE], y = y_i[id$train],
                   tau = stage1$tau, pf = d$pf, lambda = lambda,
                   lambda2 = lambda2, standardize = standardize,
                   intercept = TRUE, ...)
    pred <- as.matrix(stats::predict(fit_s, newx = d$x[id$test, , drop = FALSE], s = lambda))
    if (ncol(pred) != length(lambda)) pred <- matrix(pred, ncol = length(lambda))
    e <- y_i[id$test] - pred
    loss[s, ] <- colMeans(expectile_loss(e, stage1$tau))
    fold_rows[[s]] <- data.frame(fold = s, lambda = lambda, valid_loss = loss[s, ],
                                 n_train = length(id$train), n_valid = length(id$test),
                                 row.names = NULL)
    nested_rows[[s]] <- cbind(fold = s, nested$info)
    nested_path_rows[[s]] <- cbind(fold = s, nested$path)
  }

  cvm <- colMeans(loss, na.rm = TRUE)
  cvsd <- apply(loss, 2, stats::sd, na.rm = TRUE) / sqrt(nrow(loss))
  picked <- choose_lambda(lambda, cvm, cvsd)
  lambda_value <- picked[[lambda_choice]]
  cf <- stats::coef(full_fit, s = lambda_value)
  beta <- as.numeric(cf)
  names(beta) <- rownames(cf) %||% names(beta)
  yhat <- as.numeric(stats::predict(full_fit, newx = full_d$x, s = lambda_value))
  active <- names(beta)[names(beta) != "(Intercept)" & abs(beta) > zero_tol]
  nzero <- lasso_nzero(full_fit, lambda)
  path <- data.frame(lambda = lambda, cvm = cvm, cvsd = cvsd, nzero = nzero,
                     selected = abs(lambda - lambda_value) < 1e-12,
                     cv_method = paste0(cv, "_nested_residual"), row.names = NULL)

  structure(
    c(
      list(fit = list(ernet.fit = full_fit), beta = beta, fitted = yhat,
           resid = y_i - yhat, selected_terms = active, lambda = lambda_value,
           lambda_choice = lambda_choice, cv = cv, path = path,
           fold_loss = do.call(rbind, fold_rows),
           nested_stage1 = do.call(rbind, nested_rows),
           nested_stage1_path = do.call(rbind, nested_path_rows),
           tau = stage1$tau),
      picked
    ),
    class = "expectsee_lasso"
  )
}

fit_expectile_lasso <- function(y, x, tau, lambda_choice = c("lambda.min", "lambda.1se"),
                                penalty_factor = NULL, foldid = NULL, nfolds = 5,
                                standardize = TRUE, intercept = TRUE, lambda = NULL,
                                lambda2 = 0, zero_tol = 1e-10,
                                cv = c("iid_kfold", "fixed_window", "expanding_window"),
                                ts_splits = NULL, initial = NULL, assess = NULL,
                                step = NULL, ...) {
  lambda_choice <- match.arg(lambda_choice)
  cv <- match.arg(cv)
  check_tau(tau)
  x <- as.matrix(x)
  y <- as.numeric(y)
  if (nrow(x) != length(y)) stop("x and y have incompatible dimensions")
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))

  penalty_factor <- penalty_factor %||% rep(1, ncol(x))
  if (cv == "iid_kfold") {
    args <- list(
      x = x, y = y, tau = tau, pf = penalty_factor, lambda2 = lambda2,
      nfolds = nfolds, standardize = standardize, intercept = intercept, ...
    )
    if (!is.null(foldid)) args$foldid <- foldid
    if (!is.null(lambda)) args$lambda <- lambda
    cv_fit <- do.call(sales_call("cv.ernet"), args)
    cf <- stats::coef(cv_fit, s = lambda_choice)
    yhat <- as.numeric(stats::predict(cv_fit, newx = x, s = lambda_choice))
  } else {
    ts_splits <- ts_splits %||% see_ts_splits(
      nrow(x), initial = initial %||% floor(nrow(x) / 2),
      assess = assess %||% max(1, floor(nrow(x) / 10)),
      step = step %||% (assess %||% max(1, floor(nrow(x) / 10))),
      type = cv
    )
    cv_fit <- cv_expectile_lasso_ts(y, x, tau, ts_splits, penalty_factor,
                                    lambda, lambda2, standardize, intercept, ...)
    cf <- stats::coef(cv_fit$ernet.fit, s = cv_fit[[lambda_choice]])
    yhat <- as.numeric(stats::predict(cv_fit$ernet.fit, newx = x,
                                      s = cv_fit[[lambda_choice]]))
  }
  beta <- as.numeric(cf)
  terms <- rownames(cf) %||% names(beta)
  names(beta) <- terms

  resid <- y - yhat
  active <- terms[terms != "(Intercept)" & abs(beta) > zero_tol]
  lambda_value <- cv_fit[[lambda_choice]]
  path <- data.frame(
    lambda = cv_fit$lambda, cvm = cv_fit$cvm, cvsd = cv_fit$cvsd,
    nzero = cv_fit$nzero, selected = abs(cv_fit$lambda - lambda_value) < 1e-12,
    cv_method = cv,
    row.names = NULL
  )

  structure(
    list(
      fit = cv_fit, beta = beta, fitted = yhat, resid = resid,
      selected_terms = active, lambda = lambda_value,
      lambda_choice = lambda_choice, cv = cv, path = path, tau = tau
    ),
    class = "expectsee_lasso"
  )
}

make_active_design <- function(x_raw, selected) {
  cbind(const = 1, x_raw[, selected, drop = FALSE])
}

active_lasso_fit <- function(lasso, y, x, selected) {
  selected <- unique(selected)
  x_active <- make_active_design(x, selected)
  b0 <- if ("(Intercept)" %in% names(lasso$beta)) {
    unname(lasso$beta["(Intercept)"])
  } else if ("const" %in% names(lasso$beta)) {
    unname(lasso$beta["const"])
  } else {
    0
  }
  beta <- c(const = b0, lasso$beta[selected])
  fitted <- as.numeric(x_active %*% beta)
  lasso$beta <- beta
  lasso$fitted <- fitted
  lasso$resid <- y - fitted
  lasso$selected_terms <- selected
  list(fit = lasso, x = x_active)
}

selected_coefficient_fit <- function(lasso, y, x, selected, tau,
                                     coefficient_fit = c("lasso", "post_lasso")) {
  coefficient_fit <- match.arg(coefficient_fit)
  if (coefficient_fit == "lasso") return(active_lasso_fit(lasso, y, x, selected))

  x_active <- make_active_design(x, unique(selected))
  check_fixed_design(x_active, "post-Lasso refit")
  fit <- fit_expectile(y, x_active, tau)
  names(fit$beta) <- colnames(x_active)
  list(fit = fit, x = x_active)
}

make_lasso_coef <- function(stage, equation, beta, map, se = NA_real_,
                            coefficient_fit = "lasso") {
  terms <- names(beta)
  map_id <- match(terms, map$term)
  data.frame(
    stage = stage, coefficient_fit = coefficient_fit,
    equation = equation, term = terms,
    term_type = ifelse(terms == "const", "intercept", "lag"),
    source_equation = c(NA_character_, map$source_equation[map_id[-1]]),
    estimate = beta, sandwich_se = se,
    t = beta / se,
    p = 2 * stats::pnorm(abs(beta / se), lower.tail = FALSE),
    row.names = NULL
  )
}

hdsals <- function(y, spec = "ABSQ", tau = 0.99855, q = 1,
                   lambda_choice = c("lambda.min", "lambda.1se"),
                   foldid = NULL, nfolds = 5, standardize = TRUE,
                   own_only = FALSE,
                   coefficient_fit = c("lasso", "post_lasso"), hac_lag = 0,
                   cv = c("iid_kfold", "fixed_window", "expanding_window"),
                   ts_splits = NULL, initial = NULL, assess = NULL,
                   step = NULL, lambda = NULL, lambda2 = 0, zero_tol = 1e-10, ...) {
  lambda_choice <- match.arg(lambda_choice)
  coefficient_fit <- match.arg(coefficient_fit)
  cv <- match.arg(cv)
  design <- see_build_design(y, spec = spec, q = q, include_intercept = FALSE)
  yy <- design$y
  equations <- colnames(yy)
  k <- ncol(yy)

  coef_rows <- gof_rows <- diag_rows <- support_rows <- paths <- vector("list", k)
  fits <- selection_fits <- x_list <- vcov_list <- term_names <- term_maps <- vector("list", k)
  fitted_mat <- resid_mat <- matrix(NA_real_, nrow(yy), k, dimnames = dimnames(yy))

  for (i in seq_len(k)) {
    keep <- if (own_only) design$colmap$source_index == i else rep(TRUE, ncol(design$x))
    x_lasso <- design$x[, keep, drop = FALSE]
    map_i <- design$colmap[keep, , drop = FALSE]
    y_i <- yy[, i]

    # Stage 1 HDSALS: lasso selects lag regressors from the SEE candidate set.
    lasso <- fit_expectile_lasso(y_i, x_lasso, tau, lambda_choice, foldid = foldid,
                                 nfolds = nfolds, standardize = standardize,
                                 cv = cv, ts_splits = ts_splits, initial = initial,
                                 assess = assess, step = step, lambda = lambda,
                                 lambda2 = lambda2, zero_tol = zero_tol, ...)
    selected <- intersect(lasso$selected_terms, colnames(x_lasso))
    active <- selected_coefficient_fit(lasso, y_i, x_lasso, selected, tau,
                                       coefficient_fit)
    fit <- active$fit
    x_final <- active$x
    beta <- fit$beta
    if (coefficient_fit == "post_lasso") {
      sw <- expectile_sandwich(y_i, x_final, beta, tau, hac_lag, names(beta))
      se <- sw$se
      vc <- sw$vcov
    } else {
      se <- stats::setNames(rep(NA_real_, length(beta)), names(beta))
      vc <- matrix(NA_real_, length(beta), length(beta),
                   dimnames = list(names(beta), names(beta)))
    }

    fitted_mat[, i] <- fit$fitted
    resid_mat[, i] <- fit$resid
    fits[[i]] <- fit
    selection_fits[[i]] <- lasso
    x_list[[i]] <- x_final
    vcov_list[[i]] <- vc
    term_names[[i]] <- names(beta)
    term_maps[[i]] <- rbind(
      data.frame(term = "const", source_equation = NA_character_, source_index = NA_integer_,
                 lag = NA_integer_, transform = "intercept"),
      map_i[match(selected, map_i$term), , drop = FALSE]
    )

    coef_rows[[i]] <- make_lasso_coef(
      "HDSALS", equations[i], beta, term_maps[[i]], se, coefficient_fit
    )
    gof_rows[[i]] <- cbind(
      stage = "HDSALS", coefficient_fit = coefficient_fit,
      equation = equations[i], tau = tau, spec = design$spec,
       see_gof(y_i, fitted_mat[, i], ncol(x_final), tau,
               p_effective = qr(x_final)$rank)
    )
    diag_rows[[i]] <- cbind(stage = "HDSALS", coefficient_fit = coefficient_fit,
                            equation = equations[i], residual_diagnostics(resid_mat[, i]))
    paths[[i]] <- cbind(stage = "HDSALS", coefficient_fit = coefficient_fit,
                        equation = equations[i], lasso$path)

    support_rows[[i]] <- if (length(selected)) {
      cbind(stage = "HDSALS", coefficient_fit = coefficient_fit,
            equation = equations[i],
            map_i[match(selected, map_i$term), , drop = FALSE])
    } else {
      data.frame()
    }
  }

  structure(
    list(
      call = match.call(), tau = tau, spec = design$spec, q = q,
      own_only = own_only, coefficient_fit = coefficient_fit,
      coefficient_fit_label = if (coefficient_fit == "lasso") "Lasso" else "post-Lasso",
      inference_available = FALSE,
      coefficient_se_available = coefficient_fit == "post_lasso",
      selection_adjusted_inference = FALSE, se = if (coefficient_fit == "post_lasso") "sandwich" else "none",
      y = yy, dates = design$dates,
      design = design, x_raw = design$x, x_list = x_list, fits = fits,
      selection_fits = selection_fits,
      fitted = fitted_mat, resid = resid_mat, coef = do.call(rbind, coef_rows),
      gof = do.call(rbind, gof_rows), diagnostics = do.call(rbind, diag_rows),
      support = do.call(rbind, support_rows), lambda_path = do.call(rbind, paths),
      vcov_sandwich = vcov_list, term_names = term_names, term_maps = term_maps,
      tuning = list(
        lambda_choice = lambda_choice, nfolds = nfolds,
        standardize = standardize, cv = cv, initial = initial,
        assess = assess, step = step, lambda2 = lambda2,
        zero_tol = zero_tol
      )
    ),
    class = "expectsee_hdsals"
  )
}

hdsfgals <- function(stage1, lambda_choice = c("lambda.min", "lambda.1se"),
                     foldid = NULL, nfolds = 5, standardize = TRUE,
                     coefficient_fit = NULL, se = NULL,
                     hac_lag = 0, cv = c("iid_kfold", "fixed_window", "expanding_window"),
                     ts_splits = NULL, initial = NULL, assess = NULL, step = NULL,
                     lambda = NULL, lambda2 = 0, zero_tol = 1e-10,
                     nested_residual_cv = NULL,
                     stage2_lag_candidates = c("full", "stage1_selected"),
                     stage2_penalty = c("all", "stage1_lags_unpenalized"), ...) {
  lambda_choice <- match.arg(lambda_choice)
  if (!inherits(stage1, "expectsee_hdsals")) {
    stop("stage1 must be an HDSALS fit", call. = FALSE)
  }
  stage1_mode <- stage1$coefficient_fit %||% "lasso"
  coefficient_fit <- match.arg(coefficient_fit %||% stage1_mode,
                               c("lasso", "post_lasso"))
  if (coefficient_fit != stage1_mode) {
    stop("HDSALS and HDSFGALS must use the same coefficient_fit", call. = FALSE)
  }
  se <- if (is.null(se)) {
    if (coefficient_fit == "post_lasso") "generated" else "none"
  } else {
    match.arg(se, c("none", "generated", "sandwich"))
  }
  if (coefficient_fit == "lasso" && se != "none") {
    warning("High-dimensional Lasso inference is unavailable; `se` is ignored.", call. = FALSE)
    se <- "none"
  }
  cv <- match.arg(cv)
  stage2_lag_candidates <- match.arg(stage2_lag_candidates)
  stage2_penalty <- match.arg(stage2_penalty)
  nested_residual_cv <- nested_residual_cv %||% (cv != "iid_kfold")
  yy <- stage1$y
  e_hat <- stage1$resid
  equations <- colnames(yy)
  k <- ncol(yy)

  coef_rows <- gof_rows <- diag_rows <- edge_rows <- paths <- vector("list", k)
  stage2_support_rows <- stage2_candidate_rows <- vector("list", k)
  fold_rows <- nested_rows <- nested_path_rows <- vector("list", k)
  fits <- selection_fits <- x_list <- vcov_sand <- vcov_gen <- term_names <- term_maps <- vector("list", k)
  fitted_mat <- resid_mat <- matrix(NA_real_, nrow(yy), k, dimnames = dimnames(yy))

  for (i in seq_len(k)) {
    stage1_selected_lags <- setdiff(stage1$term_names[[i]], "const")
    lag_terms <- if (stage2_lag_candidates == "full") {
      colnames(stage1$x_raw)
    } else {
      stage1_selected_lags
    }
    forced_lags <- if (stage2_penalty == "stage1_lags_unpenalized") {
      stage1_selected_lags
    } else {
      character()
    }
    d <- stage2_design(stage1, i, e_hat, lag_terms, forced_lags)
    z_lasso <- d$x
    pf <- d$pf
    stage2_candidate_rows[[i]] <- data.frame(
      equation = equations[i],
      coefficient_fit = coefficient_fit,
      lag_candidate_rule = stage2_lag_candidates,
      penalty_rule = stage2_penalty,
      n_lag_candidates = length(d$lag_terms),
      n_forced_lags = length(d$forced_lags),
      n_residual_candidates = length(d$residual_terms),
      row.names = NULL
    )

    lasso <- if (nested_residual_cv) {
      fit_hdsfgals_lasso_nested(stage1, i, lambda_choice, standardize, lambda,
                                lambda2, zero_tol, cv, ts_splits, initial,
                                assess, step, nfolds, stage2_lag_candidates,
                                d$lag_terms, d$forced_lags, 1,
                                coefficient_fit, ...)
    } else {
      fit_expectile_lasso(yy[, i], z_lasso, stage1$tau, lambda_choice,
                          penalty_factor = pf, foldid = foldid, nfolds = nfolds,
                          standardize = standardize, cv = cv, ts_splits = ts_splits,
                          initial = initial, assess = assess, step = step,
                          lambda = lambda, lambda2 = lambda2, zero_tol = zero_tol, ...)
    }
    selected_lag_terms <- d$lag_terms[d$lag_terms %in% lasso$selected_terms]
    selected_lags_final <- d$lag_terms[d$lag_terms %in% unique(c(d$forced_lags, selected_lag_terms))]
    selected_resid_terms <- grep("^resid_", lasso$selected_terms, value = TRUE)
    selected_resid_sources <- match(sub("^resid_", "", selected_resid_terms), equations)
    resid_ok <- !is.na(selected_resid_sources)
    selected_resid_terms <- selected_resid_terms[resid_ok]
    selected_resid_sources <- selected_resid_sources[resid_ok]

    active_terms <- c(selected_lags_final, selected_resid_terms)
    active <- selected_coefficient_fit(
      lasso, yy[, i], z_lasso, active_terms, stage1$tau, coefficient_fit
    )
    fit <- active$fit
    z_final <- active$x
    beta <- fit$beta
    if (coefficient_fit == "post_lasso") {
      sw <- expectile_sandwich(yy[, i], z_final, beta, stage1$tau,
                               hac_lag, names(beta))
      z_lag_final <- make_active_design(stage1$x_raw, selected_lags_final)
      stage2_x_list <- stage1$x_list
      stage2_x_list[[i]] <- z_lag_final
      gen <- generated_se_i(
        i, stage2_x_list, e_hat, fit$resid, beta, stage1$tau,
        hac_lag, names(beta), residual_sources = selected_resid_sources,
        first_stage_x = stage1$x_list
      )
      sw_se <- sw$se; gen_se <- gen$se
      vc_s <- sw$vcov; vc_g <- gen$vcov
    } else {
      sw_se <- gen_se <- stats::setNames(rep(NA_real_, length(beta)), names(beta))
      vc_s <- vc_g <- matrix(NA_real_, length(beta), length(beta),
                             dimnames = list(names(beta), names(beta)))
    }
    se_use <- switch(se, sandwich = sw_se, generated = gen_se,
                     none = stats::setNames(rep(NA_real_, length(beta)), names(beta)))

    fitted_mat[, i] <- fit$fitted
    resid_mat[, i] <- fit$resid
    fits[[i]] <- fit
    selection_fits[[i]] <- lasso
    x_list[[i]] <- z_final
    vcov_sand[[i]] <- vc_s
    vcov_gen[[i]] <- vc_g
    term_names[[i]] <- names(beta)

    lag_map <- stage1$design$colmap[match(selected_lags_final, stage1$design$colmap$term), , drop = FALSE]
    term_type <- c("intercept", lag_map$transform, rep("generated_residual", length(selected_resid_terms)))
    source_eq <- c(NA_character_, lag_map$source_equation, equations[selected_resid_sources])
    term_maps[[i]] <- data.frame(
      term = names(beta), source_equation = source_eq,
      source_index = c(NA_integer_, lag_map$source_index, selected_resid_sources),
      lag = c(NA_integer_, lag_map$lag, rep(NA_integer_, length(selected_resid_terms))),
      transform = term_type, stringsAsFactors = FALSE
    )
    coef_rows[[i]] <- data.frame(
      stage = "HDSFGALS", coefficient_fit = coefficient_fit,
      equation = equations[i], term = names(beta),
      term_type = ifelse(term_type == "intercept", "intercept",
                         ifelse(term_type == "generated_residual", "generated_residual", "lag")),
      source_equation = source_eq, estimate = beta,
      sandwich_se = sw_se, generated_se = gen_se,
      t = beta / se_use,
      p = 2 * stats::pnorm(abs(beta / se_use), lower.tail = FALSE),
      row.names = NULL
    )

    gof_rows[[i]] <- cbind(
      stage = "HDSFGALS", coefficient_fit = coefficient_fit,
      equation = equations[i], tau = stage1$tau, spec = stage1$spec,
       see_gof(yy[, i], fitted_mat[, i], ncol(z_final), stage1$tau,
               p_effective = qr(z_final)$rank)
    )
    diag_rows[[i]] <- cbind(stage = "HDSFGALS", coefficient_fit = coefficient_fit,
                            equation = equations[i], residual_diagnostics(resid_mat[, i]))
    paths[[i]] <- cbind(stage = "HDSFGALS", coefficient_fit = coefficient_fit,
                        equation = equations[i], lasso$path)
    fold_rows[[i]] <- if (!is.null(lasso$fold_loss)) cbind(equation = equations[i], lasso$fold_loss) else data.frame()
    nested_rows[[i]] <- if (!is.null(lasso$nested_stage1)) cbind(target_equation = equations[i], lasso$nested_stage1) else data.frame()
    nested_path_rows[[i]] <- if (!is.null(lasso$nested_stage1_path)) {
      cbind(target_equation = equations[i], lasso$nested_stage1_path)
    } else {
      data.frame()
    }
    lag_support <- if (length(selected_lags_final)) {
      data.frame(
        stage = "HDSFGALS", coefficient_fit = coefficient_fit,
        equation = equations[i],
        component = "lag",
        selected_by = ifelse(selected_lags_final %in% d$forced_lags, "forced", "lasso"),
        lag_map,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame()
    }
    residual_support <- if (length(selected_resid_sources)) {
      data.frame(
        stage = "HDSFGALS", coefficient_fit = coefficient_fit,
        equation = equations[i],
        component = "generated_residual", selected_by = "lasso",
        term = selected_resid_terms,
        source_equation = equations[selected_resid_sources],
        source_index = selected_resid_sources,
        lag = NA_integer_, transform = "generated_residual",
        row.names = NULL
      )
    } else {
      data.frame()
    }
    stage2_support_rows[[i]] <- rbind(lag_support, residual_support)
    edge_rows[[i]] <- if (length(selected_resid_sources)) {
      data.frame(
        coefficient_fit = coefficient_fit,
        source = equations[selected_resid_sources], target = equations[i],
        term = selected_resid_terms, estimate = beta[selected_resid_terms],
        row.names = NULL
      )
    } else {
      data.frame()
    }
  }

  structure(
    list(
      call = match.call(), tau = stage1$tau, spec = stage1$spec, q = stage1$q,
      first_stage = stage1, coefficient_fit = coefficient_fit,
      coefficient_fit_label = if (coefficient_fit == "lasso") "Lasso" else "post-Lasso",
      inference_available = FALSE,
      coefficient_se_available = coefficient_fit == "post_lasso",
      selection_adjusted_inference = FALSE, y = yy, dates = stage1$dates,
      fits = fits, selection_fits = selection_fits,
      x_list = x_list, fitted = fitted_mat, resid = resid_mat,
      coef = do.call(rbind, coef_rows), gof = do.call(rbind, gof_rows),
      diagnostics = do.call(rbind, diag_rows), support = stage1$support,
      stage2_support = do.call(rbind, stage2_support_rows),
      stage2_candidates = do.call(rbind, stage2_candidate_rows),
      edges = do.call(rbind, edge_rows), lambda_path = do.call(rbind, paths),
      cv_fold_loss = do.call(rbind, fold_rows),
      nested_stage1 = do.call(rbind, nested_rows),
      nested_stage1_path = do.call(rbind, nested_path_rows),
      vcov_sandwich = vcov_sand, vcov_generated = vcov_gen,
      term_names = term_names, term_maps = term_maps, se = se,
      nested_residual_cv = nested_residual_cv,
      stage2_lag_candidates = stage2_lag_candidates,
      stage2_penalty = stage2_penalty
    ),
    class = "expectsee_hdsfgals"
  )
}

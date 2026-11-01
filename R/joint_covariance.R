joint_fixed_stage <- function(obj, stage) {
  if (inherits(obj, c("expectsee_hdsals", "expectsee_hdsfgals")) ||
      any(obj$coefficient_fit %in% c("lasso", "post_lasso"))) {
    stop("Joint covariance is unavailable for high-dimensional or post-Lasso fits.",
         call. = FALSE)
  }
  expected_class <- if (stage == "SALS") "expectsee_sals" else "expectsee_sfgals"
  if (!inherits(obj, expected_class)) {
    stop("fit must contain a fixed-dimensional ", stage, " object", call. = FALSE)
  }
  if (!is.numeric(obj$tau) || length(obj$tau) != 1L ||
      !is.finite(obj$tau) || obj$tau <= 0 || obj$tau >= 1) {
    stop("fit has an invalid tau", call. = FALSE)
  }
  y <- obj$y
  e <- obj$resid
  if (!is.matrix(y) || !is.numeric(y) || any(!is.finite(y)) ||
      nrow(y) < 2L || ncol(y) < 1L ||
      !is.matrix(e) || !is.numeric(e) || any(!is.finite(e)) ||
      !identical(dim(y), dim(e)) || !identical(dimnames(y), dimnames(e))) {
    stop(stage, " needs finite, aligned outcome and residual matrices", call. = FALSE)
  }
  equations <- colnames(y)
  if (is.null(equations) || anyNA(equations) || any(!nzchar(equations)) ||
      anyDuplicated(equations)) {
    stop("fit needs unique equation names", call. = FALSE)
  }
  k <- ncol(y)
  if (!is.list(obj$x_list) || length(obj$x_list) != k ||
      !is.list(obj$term_names) || length(obj$term_names) != k ||
      !is.data.frame(obj$coef) ||
      !all(c("equation", "term", "estimate") %in% names(obj$coef)) ||
      anyNA(obj$coef$equation) || anyNA(obj$coef$term) ||
      !is.numeric(obj$coef$estimate) || any(!is.finite(obj$coef$estimate)) ||
      any(!obj$coef$equation %in% equations)) {
    stop(stage, " has malformed design or coefficient metadata", call. = FALSE)
  }
  if (!is.null(obj$coef$converged) &&
      (anyNA(obj$coef$converged) || !all(obj$coef$converged))) {
    stop(stage, " has unconverged coefficient estimates", call. = FALSE)
  }
  beta <- vector("list", k)
  index <- vector("list", k)
  for (i in seq_len(k)) {
    x <- obj$x_list[[i]]
    terms <- obj$term_names[[i]]
    label <- paste(stage, "equation", equations[i])
    if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) ||
        nrow(x) != nrow(y) || ncol(x) < 1L || ncol(x) >= nrow(x) ||
        !is.character(terms) || anyNA(terms) || any(!nzchar(terms)) ||
        anyDuplicated(terms) || !identical(colnames(x), terms)) {
      stop(label, " has an invalid or misaligned design", call. = FALSE)
    }
    h <- crossprod(x, x * omega_tau(e[, i], obj$tau)) / nrow(x)
    if (qr(x)$rank != ncol(x) || qr(h)$rank != ncol(x) ||
        is.null(tryCatch(solve(h), error = function(e) NULL))) {
      stop(label, " needs a full-rank design and nonsingular Hessian", call. = FALSE)
    }
    rows <- obj$coef[obj$coef$equation == equations[i], , drop = FALSE]
    if (nrow(rows) != length(terms) || anyDuplicated(rows$term) ||
        !setequal(rows$term, terms)) {
      stop(label, " has inconsistent coefficient terms", call. = FALSE)
    }
    beta[[i]] <- rows$estimate[match(terms, rows$term)]
    implied <- as.numeric(y[, i] - x %*% beta[[i]])
    if (max(abs(implied - e[, i])) > 1e-8 * max(1, abs(y[, i]), abs(e[, i]))) {
      stop(label, " has inconsistent coefficients, residuals, or sample order", call. = FALSE)
    }
    index[[i]] <- data.frame(stage = stage, equation = equations[i], term = terms,
                            stringsAsFactors = FALSE)
  }
  list(beta = beta, index = do.call(rbind, index))
}

see_vcov_joint <- function(fit, hac_lag = 0) {
  obj <- active_stage(fit)
  if (inherits(obj, c("expectsee_hdsals", "expectsee_hdsfgals")) ||
      any(obj$coefficient_fit %in% c("lasso", "post_lasso"))) {
    stop("Joint covariance is unavailable for high-dimensional or post-Lasso fits.",
         call. = FALSE)
  }
  has_stage2 <- inherits(obj, "expectsee_sfgals")
  first <- if (has_stage2) obj$first_stage else obj
  info1 <- joint_fixed_stage(first, "SALS")
  n <- nrow(first$y)
  if (!is.numeric(hac_lag) || length(hac_lag) != 1L || !is.finite(hac_lag) ||
      hac_lag < 0 || hac_lag != floor(hac_lag) || hac_lag >= n) {
    stop("hac_lag must be an integer between 0 and n - 1", call. = FALSE)
  }
  # Full-rank validation above makes the legacy influence helper's inverse exact.
  influence <- do.call(cbind, stage1_influence(first$x_list, first$resid, first$tau))
  index <- info1$index
  coefficients <- unlist(info1$beta, use.names = FALSE)
  if (has_stage2) {
    info2 <- joint_fixed_stage(obj, "SFGALS")
    if (!identical(obj$y, first$y) || !identical(obj$dates, first$dates) ||
        !identical(obj$tau, first$tau)) {
      stop("SALS and SFGALS must use the same outcomes, sample order, dates, and tau",
           call. = FALSE)
    }
    k <- ncol(first$y)
    influence2 <- vector("list", k)
    for (i in seq_len(k)) {
      other <- setdiff(seq_len(k), i)
      z <- cbind(first$x_list[[i]], first$resid[, other, drop = FALSE])
      terms <- c(first$term_names[[i]], paste0("resid_", colnames(first$y)[other]))
      if (!identical(obj$term_names[[i]], terms) ||
          !isTRUE(all.equal(unname(obj$x_list[[i]]), unname(z), tolerance = 0))) {
        stop("SFGALS design must contain the aligned SALS design and generated residuals",
             call. = FALSE)
      }
      gen <- generated_se_i(i, first$x_list, first$resid, obj$resid[, i],
                            info2$beta[[i]], first$tau, hac_lag = 0)
      # generated_se_i returns a corrected score; transform it to coefficient IF.
      influence2[[i]] <- gen$influence %*% gen$h_inv
    }
    influence <- cbind(influence, do.call(cbind, influence2))
    index <- rbind(index, info2$index)
    coefficients <- c(coefficients, unlist(info2$beta, use.names = FALSE))
  }
  if (any(!is.finite(influence))) stop("Coefficient influences are not finite", call. = FALSE)
  labels <- make.unique(paste(index$stage, index$equation, index$term, sep = "::"))
  colnames(influence) <- names(coefficients) <- labels
  vcov_sqrt_n <- sym(if (hac_lag > 0) meat_hac(influence, hac_lag) else meat_iid(influence))
  list(vcov = vcov_sqrt_n / n, vcov_sqrt_n = vcov_sqrt_n,
       vcov_sqrt_t = vcov_sqrt_n, influence = influence, index = index,
       coefficients = coefficients, n = n, hac_lag = as.integer(hac_lag))
}

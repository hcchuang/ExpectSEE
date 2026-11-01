#' Impulse responses for expectiles (IRFE)
#'
#' A one-time Stage-2 innovation d moves all equations on impact by A d.
#' center = "zero" starts both q-lag histories at zero and retains the original G:
#' b_1 = G(0), x_1 = G(0) + A d; thereafter b and x both iterate G, and IRFE = x - b.
#' center = "fixed_point" starts both histories at a locally stable root ybar = G(ybar).
#' Future innovations are zero. No intercept is removed in either mode.
#' Fixed-dimensional intervals use the joint SALS-SFGALS covariance and propagate
#' derivatives through both complete paths. HDSFGALS supports point estimates only.
see_irfe <- function(fit, source = NULL, shock = 0.01, horizon = 12L,
                     sign = c("both", "positive", "negative"), ci = FALSE, level = 0.95,
                     hac_lag = 4L, rcond_tol = 1e-8, divergence = c("error", "record"),
                     fd_step = 1e-6, center = c("zero", "fixed_point"), fixed_point = NULL,
                     fp_tol = 1e-12, ci_gradient = c("recursive", "finite_difference")) {
  sign <- match.arg(sign)
  divergence <- match.arg(divergence)
  center <- match.arg(center)
  ci_gradient <- match.arg(ci_gradient)
  if (length(fp_tol) != 1L || !is.numeric(fp_tol) || !is.finite(fp_tol) || fp_tol <= 0) {
    stop("fp_tol must be one positive finite number", call. = FALSE)
  }
  if (center == "zero" && !is.null(fixed_point)) {
    stop("fixed_point is used only with center = \"fixed_point\"", call. = FALSE)
  }
  if (length(shock) != 1L || !is.numeric(shock) || !is.finite(shock) || shock <= 0) {
    stop("shock must be one positive finite number; use sign for negative innovations", call. = FALSE)
  }
  if (length(horizon) != 1L || !is.numeric(horizon) || !is.finite(horizon) || horizon < 1 ||
      horizon != as.integer(horizon)) {
    stop("horizon must be a positive integer (h = 1 is the impact)", call. = FALSE)
  }
  if (!is.logical(ci) || length(ci) != 1L || is.na(ci)) stop("ci must be TRUE or FALSE", call. = FALSE)
  if (length(level) != 1L || !is.numeric(level) || !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be in (0, 1)", call. = FALSE)
  }
  if (length(fd_step) != 1L || !is.numeric(fd_step) || !is.finite(fd_step) || fd_step <= 0) {
    stop("fd_step must be one positive finite number", call. = FALSE)
  }
  horizon <- as.integer(horizon)
  system <- see_system_map(fit, rcond_tol)
  eq <- system$equations
  k <- length(eq)
  sources <- if (is.null(source)) eq else if (is.numeric(source)) eq[source] else source
  if (!length(sources) || anyNA(sources) || !all(sources %in% eq) || anyDuplicated(sources)) {
    stop("source must name or index system equations", call. = FALSE)
  }
  signs <- switch(sign, both = c(1, -1), positive = 1, negative = -1)
  cases <- expand.grid(sign = signs, source = sources, stringsAsFactors = FALSE)[c("source", "sign")]
  rest <- .irfe_rest(system, center, fixed_point, fp_tol)
  paths <- lapply(seq_len(nrow(cases)), function(r)
    .irfe_path(system, cases$source[r], cases$sign[r] * shock, horizon, divergence, rest$base, center == "fixed_point"))
  g0 <- .system_next(system, matrix(0, system$q, k, dimnames = list(NULL, eq)), rep(0, k))
  out <- do.call(rbind, lapply(seq_len(nrow(cases)), function(r) {
    p <- paths[[r]]
    data.frame(source = cases$source[r], sign = cases$sign[r], shock = cases$sign[r] * shock,
               target = rep(eq, horizon), horizon = rep(seq_len(horizon), each = k),
               irfe = as.vector(t(p$path)),
               baseline = as.vector(t(p$baseline)), shocked = as.vector(t(p$shocked)),
               initial_transition = rep(unname(g0), horizon), removed_intercept = NA_real_,
               path_boundary = rep(p$boundary, each = k),
               path_status = p$status, first_nonfinite_horizon = p$first_bad,
               divergent_path = p$divergent_path,
               stringsAsFactors = FALSE)
  }))
  se <- if (ci) {
    if (ci_gradient == "recursive")
      .irfe_se_recursive(fit, system, cases, shock, horizon, hac_lag, center, rest$base)
    else .irfe_se(fit, system, cases, shock, horizon, hac_lag, fd_step, center, rest$base,
                  fp_tol, rcond_tol)
  } else NA_real_
  z <- stats::qnorm(1 - (1 - level) / 2)
  meta <- data.frame(estimator = system$estimator, coefficient_fit = system$coefficient_fit,
                     tau = system$tau, spec = system$spec, q = system$q, stringsAsFactors = FALSE)
  out <- cbind(meta[rep(1L, nrow(out)), , drop = FALSE], out, row.names = NULL)
  out$se <- se
  out$ci_status <- if (!ci) "not_requested" else ifelse(!is.finite(out$irfe), "path_nonfinite",
    ifelse(out$path_boundary, "switching_boundary", ifelse(!is.finite(out$se), "derivative_nonfinite",
      ifelse(out$se <= 0, "degenerate_variance", "available"))))
  out$se[out$ci_status == "switching_boundary"] <- NA_real_
  out$lower <- out$irfe - z * out$se
  out$upper <- out$irfe + z * out$se
  out$p_value <- 2 * stats::pnorm(-abs(out$irfe / out$se))
  out[out$ci_status != "available", c("lower", "upper", "p_value")] <- NA_real_
  out$level <- if (ci) level else NA_real_
  out$hac_lag <- if (ci) as.integer(hac_lag) else NA_integer_
  out$ci_method <- if (ci) ifelse(out$horizon == 1L, "analytic_dA",
                                  if (ci_gradient == "recursive") "recursive_paired_path" else "central_difference") else NA_character_
  out$rcond <- system$diagnostics$rcond
  out$inverse_norm_inf <- system$diagnostics$inverse_norm_inf
  out$response_type <- if (center == "zero") "irfe_zero_initial_state" else "irfe_fixed_point_centered"
  out$inference_available <- out$ci_status == "available"
  out$center <- center
  out$fd_step <- if (ci && ci_gradient == "finite_difference") fd_step else NA_real_
  out$fixed_point_residual <- rest$residual
  out$fixed_point_spectral_radius <- rest$spectral_radius
  class(out) <- c("expectsee_irfe", "data.frame")
  out
}

# Evolve both full q-lag histories using the original G. The innovation is added
# only after the first transition; all later innovations are zero.
.irfe_path <- function(system, source, d, horizon, divergence = "error", base = NULL,
                       estimated_base = FALSE) {
  eq <- system$equations
  k <- length(eq)
  q <- system$q
  if (is.null(base)) base <- rep(0, k)
  bstate <- xstate <- matrix(base, q, k, byrow = TRUE, dimnames = list(NULL, eq))
  baseline <- shocked <- path <- matrix(NA_real_, horizon, k, dimnames = list(NULL, eq))
  first_bad <- NA_integer_
  divergent_path <- NA_character_
  boundary <- rep(FALSE, horizon)
  # Only active nonsquared positive/negative terms can have a first-derivative kink.
  terms <- which(system$basis$transform %in% c("plus", "minus") &
                   colSums(abs(system$reduced_coefficients)) > 0)
  at_boundary <- function(state, h) any(vapply(terms, function(i) {
    term <- system$basis[i, ]
    if (!estimated_base && term$lag >= h) return(FALSE) # fixed initial zero rows
    abs(state[q - term$lag + 1L, term$source_equation]) <= 1e-8
  }, logical(1)))
  advance <- function(state) tryCatch(.system_next(system, state, rep(0, k)), error = identity)
  for (h in seq_len(horizon)) {
    if (h > 1L) boundary[h] <- boundary[h - 1L] || at_boundary(bstate, h) || at_boundary(xstate, h)
    b <- advance(bstate)
    x <- advance(xstate)
    if (!inherits(x, "error") && h == 1L) x <- x + system$A[, source] * d
    bad_b <- inherits(b, "error") || any(!is.finite(b))
    bad_x <- inherits(x, "error") || any(!is.finite(x))
    if (bad_b || bad_x || any(!is.finite(x - b))) {
      divergent_path <- if (bad_b && bad_x) "both" else if (bad_b) "baseline" else if (bad_x) "shocked" else "response"
      if (divergence == "error") stop("IRFE ", divergent_path, " path is not finite at horizon ", h, call. = FALSE)
      first_bad <- h
      break
    }
    baseline[h, ] <- b
    shocked[h, ] <- x
    # Retain the exact analytic impact instead of subtracting nearly equal levels.
    path[h, ] <- if (h == 1L) system$A[, source] * d else x - b
    bstate <- utils::tail(rbind(bstate, b), q)
    xstate <- utils::tail(rbind(xstate, x), q)
  }
  list(path = path, baseline = baseline, shocked = shocked, boundary = boundary,
       status = if (is.na(first_bad)) "finite" else "diverged", first_bad = first_bad,
       divergent_path = divergent_path)
}

# The reduced transition at a constant history (all q lags equal to y) and its Jacobians with
# respect to each lag row of the state (rows ordered oldest to newest).
.irfe_map <- function(system) {
  eq <- system$equations
  k <- length(eq)
  q <- system$q
  rows <- function(y) matrix(y, q, k, byrow = TRUE, dimnames = list(NULL, eq))
  step <- function(state) .system_next(system, state, rep(0, k))
  lag_jac <- function(y, h = 1e-7) lapply(seq_len(q), function(r) vapply(seq_len(k), function(j) {
    up <- rows(y)
    dn <- up
    up[r, j] <- up[r, j] + h
    dn[r, j] <- dn[r, j] - h
    unname(step(up) - step(dn)) / (2 * h)
  }, numeric(k)))
  list(G = function(y) unname(step(rows(y))), lag_jac = lag_jac, k = k, q = q)
}

# Fixed point y = G(y) by damped Newton from start (one polishing step after convergence);
# NULL if Newton fails.
.irfe_newton <- function(map, start, tol, max_iter = 100L) {
  resid <- function(y) {
    f <- tryCatch(map$G(y) - y, error = function(e) NULL)
    if (is.null(f) || !all(is.finite(f))) NULL else f
  }
  newton_step <- function(y, f) {
    j <- Reduce(`+`, map$lag_jac(y)) - diag(map$k)
    tryCatch(solve(j, f), error = function(e) NULL)
  }
  y <- as.numeric(start)
  for (it in seq_len(max_iter)) {
    f <- resid(y)
    if (is.null(f)) return(NULL)
    if (max(abs(f)) < tol) {
      s <- newton_step(y, f)
      if (!is.null(s)) {
        f1 <- resid(y - s)
        if (!is.null(f1) && max(abs(f1)) <= max(abs(f))) y <- y - s
      }
      return(y)
    }
    s <- newton_step(y, f)
    if (is.null(s)) return(NULL)
    t <- 1
    repeat {
      f1 <- resid(y - t * s)
      if (!is.null(f1) && max(abs(f1)) < max(abs(f))) break
      t <- t / 2
      if (t < 1e-10) return(NULL)
    }
    y <- y - t * s
  }
  NULL
}

# Spectral radius of the companion Jacobian of the state map at the constant history y
# (below one: the fixed point is locally stable).
.irfe_spectral_radius <- function(map, y) {
  k <- map$k
  q <- map$q
  jl <- map$lag_jac(y)
  comp <- matrix(0, q * k, q * k)
  if (q > 1L) for (r in seq_len(q - 1L)) comp[(r - 1L) * k + seq_len(k), r * k + seq_len(k)] <- diag(k)
  for (r in seq_len(q)) comp[(q - 1L) * k + seq_len(k), (r - 1L) * k + seq_len(k)] <- jl[[r]]
  max(Mod(eigen(comp, only.values = TRUE)$values))
}

# Initial state: zero, or the unique locally stable fixed point of G found by Newton
# from G(0) and from the iterated map, or a supplied fixed point (verified).
.irfe_rest <- function(system, center, fixed_point, tol) {
  k <- length(system$equations)
  if (center == "zero") return(list(base = rep(0, k), residual = NA_real_, spectral_radius = NA_real_))
  map <- .irfe_map(system)
  if (!is.null(fixed_point)) {
    start <- unname(.system_vector(fixed_point, system$equations, "fixed_point"))
    r <- .irfe_newton(map, start, tol)
    if (is.null(r) || max(abs(r - start)) > 1e-6) {
      stop("fixed_point is not a fixed point of the reduced transition", call. = FALSE)
    }
    roots <- list(r)
  } else {
    starts <- list()
    g0 <- tryCatch(map$G(rep(0, k)), error = function(e) NULL)
    if (!is.null(g0) && all(is.finite(g0))) starts[[1L]] <- g0
    y <- rep(0, k)
    for (it in seq_len(5000L)) {
      y1 <- tryCatch(map$G(y), error = function(e) NULL)
      if (is.null(y1) || !all(is.finite(y1)) || max(abs(y1)) > 1e6) break
      if (max(abs(y1 - y)) < 1e-10) {
        starts[[length(starts) + 1L]] <- y1
        break
      }
      y <- y1
    }
    roots <- Filter(Negate(is.null), lapply(starts, function(s) .irfe_newton(map, s, tol)))
  }
  rho <- vapply(roots, function(r) .irfe_spectral_radius(map, r), numeric(1))
  keep <- roots[rho < 1]
  if (!length(keep)) {
    stop(if (length(roots)) sprintf("the fixed point found is not locally stable (spectral radius %.3g); ", min(rho))
         else "no fixed point of the reduced transition was found; ",
         "use center = \"zero\"", if (is.null(fixed_point)) " or supply fixed_point" else "",
         call. = FALSE)
  }
  distinct <- keep[!duplicated(lapply(keep, function(r) round(r, 8)))]
  if (length(distinct) > 1L) stop("several locally stable fixed points were found; supply fixed_point", call. = FALSE)
  r <- distinct[[1L]]
  list(base = r, residual = max(abs(map$G(r) - r)), spectral_radius = .irfe_spectral_radius(map, r))
}

# IRFE paths after replacing all SALS and SFGALS coefficients by theta (ordered as the
# see_vcov_joint() index). With center = "fixed_point", the fixed point is re-solved by Newton
# from base for the perturbed coefficients (NA paths if that fails).
.irfe_paths_theta <- function(obj, index, theta, cases, shock, horizon, center = "zero",
                              base = NULL, tol = 1e-12, rcond_tol = 1e-8) {
  set <- function(tab, rows) {
    key <- match(paste(tab$equation, tab$term, sep = "\r"),
                 paste(index$equation[rows], index$term[rows], sep = "\r"))
    if (anyNA(key)) stop("coefficient table does not match the joint covariance index", call. = FALSE)
    tab$estimate <- theta[rows][key]
    tab
  }
  stage1 <- which(index$stage == "SALS")
  obj$first_stage$coef <- set(obj$first_stage$coef, stage1)
  obj$coef <- set(obj$coef, setdiff(seq_along(theta), stage1))
  k <- ncol(obj$y)
  unavailable <- function() lapply(seq_len(nrow(cases)), function(r) matrix(NA_real_, horizon, k))
  system <- tryCatch(see_system_map(obj, rcond_tol), error = function(e) NULL)
  if (is.null(system)) return(unavailable())
  if (is.null(base)) base <- rep(0, k)
  if (center == "fixed_point") {
    base <- tryCatch(.irfe_newton(.irfe_map(system), base, tol), error = function(e) NULL)
    if (is.null(base)) return(unavailable())
  }
  lapply(seq_len(nrow(cases)), function(r)
    .irfe_path(system, cases$source[r], cases$sign[r] * shock, horizon, "record", base,
               center == "fixed_point")$path)
}

# Derivative of a basis term away from an active plus/minus switching point.
.irfe_feature_derivative <- function(x, transform) {
  switch(transform, intercept = 0, raw = 1,
         plus = as.numeric(x > 0), minus = -as.numeric(x < 0),
         plus2 = 2 * max(x, 0), minus2 = -2 * max(-x, 0))
}

# Analytic paired-path delta method for every fixed-dimensional specification.
# The joint covariance index supplies both stages and all cross-equation blocks.
.irfe_se_recursive <- function(fit, system, cases, shock, horizon, hac_lag,
                               center = "zero", base = NULL) {
  obj <- active_stage(fit)
  if (!inherits(obj, "expectsee_sfgals") ||
      any(obj$coefficient_fit %in% c("lasso", "post_lasso"))) {
    stop("IRFE intervals are available only for fixed-dimensional SFGALS fits; ",
         "high-dimensional (LASSO) IRFEs are point estimates", call. = FALSE)
  }
  v <- see_vcov_joint(obj, hac_lag = hac_lag)
  eq <- system$equations
  basis <- system$basis
  k <- length(eq)
  q <- system$q
  m <- nrow(basis)
  p <- length(v$coefficients)
  a <- system$A
  d <- system$D
  c1 <- system$stage1_coefficients
  rcoef <- system$reduced_coefficients
  dr <- array(0, c(k, m, p))
  impact <- array(0, c(k, k, p))
  ad <- a %*% d
  for (t in seq_len(p)) {
    i <- match(v$index$equation[t], eq)
    term <- v$index$term[t]
    if (v$index$stage[t] == "SALS") {
      j <- match(term, basis$term)
      if (is.na(j)) stop("joint covariance does not match the Stage-1 basis", call. = FALSE)
      dr[, j, t] <- -ad[, i]
    } else if (startsWith(term, "resid_")) {
      j <- match(sub("^resid_", "", term), eq)
      if (is.na(j) || i == j || abs(v$coefficients[t] - d[i, j]) > 1e-10)
        stop("joint covariance does not match the generated-residual coefficients", call. = FALSE)
      dr[, , t] <- tcrossprod(a[, i], rcoef[j, ] - c1[j, ])
      impact[, , t] <- tcrossprod(a[, i], a[j, ])
    } else {
      j <- match(term, basis$term)
      if (is.na(j)) stop("joint covariance does not match the Stage-2 basis", call. = FALSE)
      dr[, j, t] <- a[, i]
    }
  }
  features <- function(state) .system_features(system, state)
  direct <- function(phi) {
    out <- matrix(0, k, p)
    for (t in seq_len(p)) out[, t] <- matrix(dr[, , t], k, m) %*% phi
    out
  }
  advance <- function(state, state_jac) {
    phi <- features(state)
    jac <- direct(phi)
    for (j in seq_len(m)) {
      b <- basis[j, ]
      if (b$transform == "intercept") next
      row <- q - b$lag + 1L
      src <- b$source_index
      slope <- .irfe_feature_derivative(state[row, src], b$transform)
      if (slope != 0) jac <- jac + tcrossprod(rcoef[, j] * slope, state_jac[row, src, ])
    }
    list(value = as.vector(rcoef %*% phi), jac = jac)
  }
  if (is.null(base)) base <- rep(0, k)
  base_jac <- matrix(0, k, p)
  if (center == "fixed_point") {
    steady <- matrix(base, q, k, byrow = TRUE, dimnames = list(NULL, eq))
    state_jacobian <- matrix(0, k, k)
    for (j in seq_len(m)) {
      b <- basis[j, ]
      if (b$transform == "intercept") next
      slope <- .irfe_feature_derivative(base[b$source_index], b$transform)
      state_jacobian[, b$source_index] <-
        state_jacobian[, b$source_index] + rcoef[, j] * slope
    }
    base_jac <- tryCatch(solve(diag(k) - state_jacobian, direct(features(steady))),
                         error = function(e) matrix(NA_real_, k, p))
  }
  shift <- function(jac, next_jac) {
    out <- array(0, c(q, k, p))
    if (q > 1L) out[seq_len(q - 1L), , ] <- jac[2:q, , ]
    out[q, , ] <- next_jac
    out
  }
  se <- lapply(seq_len(nrow(cases)), function(case) {
    bstate <- xstate <- matrix(base, q, k, byrow = TRUE, dimnames = list(NULL, eq))
    bjac <- xjac <- array(rep(base_jac, each = q), c(q, k, p))
    ans <- matrix(NA_real_, horizon, k)
    for (h in seq_len(horizon)) {
      b <- advance(bstate, bjac)
      x <- advance(xstate, xjac)
      if (h == 1L) {
        source <- match(cases$source[case], eq)
        size <- cases$sign[case] * shock
        x$value <- x$value + a[, source] * size
        x$jac <- x$jac + matrix(impact[, source, ], k, p) * size
      }
      if (any(!is.finite(c(b$value, x$value, b$jac, x$jac)))) break
      gradient <- x$jac - b$jac
      ans[h, ] <- sqrt(pmax(rowSums((gradient %*% v$vcov) * gradient), 0))
      bstate <- utils::tail(rbind(bstate, b$value), q)
      xstate <- utils::tail(rbind(xstate, x$value), q)
      bjac <- shift(bjac, b$jac)
      xjac <- shift(xjac, x$jac)
    }
    ans
  })
  unlist(lapply(se, function(x) as.vector(t(x))), use.names = FALSE)
}

# Delta-method standard errors in the row order of see_irfe(). The entire paired
# path is re-evaluated at each theta, including the parameter-dependent baseline.
.irfe_se <- function(fit, system, cases, shock, horizon, hac_lag, fd_step, center = "zero",
                     base = NULL, tol = 1e-12, rcond_tol = 1e-8) {
  obj <- active_stage(fit)
  if (!inherits(obj, "expectsee_sfgals") || any(obj$coefficient_fit %in% c("lasso", "post_lasso"))) {
    stop("IRFE intervals are available only for fixed-dimensional SFGALS fits; ",
         "high-dimensional (LASSO) IRFEs are point estimates", call. = FALSE)
  }
  v <- see_vcov_joint(obj, hac_lag = hac_lag)
  eq <- system$equations
  k <- length(eq)
  idx <- which(v$index$stage == "SFGALS" & grepl("^resid_", v$index$term))
  a_eq <- v$index$equation[idx]
  b_eq <- sub("^resid_", "", v$index$term[idx])
  if (!length(idx) || max(abs(v$coefficients[idx] - system$D[cbind(a_eq, b_eq)])) > 1e-10) {
    stop("joint covariance does not match the generated-residual coefficients", call. = FALSE)
  }
  vd <- v$vcov[idx, idx, drop = FALSE]
  a <- system$A
  se <- lapply(seq_len(nrow(cases)), function(r) {
    d <- cases$sign[r] * shock
    j <- cases$source[r]
    out <- matrix(NA_real_, horizon, k)
    # h = 1: d A_ij / d D_ab = A_ia A_bj (dA = A dD A).
    out[1L, ] <- vapply(eq, function(i) {
      g <- d * a[i, a_eq] * a[b_eq, j]
      sqrt(max(drop(crossprod(g, vd %*% g)), 0))
    }, numeric(1))
    out
  })
  if (horizon > 1L) {
    theta <- v$coefficients
    n_par <- length(theta)
    jac <- lapply(seq_len(nrow(cases)), function(r) array(0, c(horizon - 1L, k, n_par)))
    for (p in seq_len(n_par)) {
      step <- fd_step * max(1, abs(theta[p]))
      up <- theta; up[p] <- theta[p] + step
      down <- theta; down[p] <- theta[p] - step
      plus <- .irfe_paths_theta(obj, v$index, up, cases, shock, horizon, center, base, tol, rcond_tol)
      minus <- .irfe_paths_theta(obj, v$index, down, cases, shock, horizon, center, base, tol, rcond_tol)
      for (r in seq_len(nrow(cases))) {
        jac[[r]][, , p] <- (plus[[r]][-1L, , drop = FALSE] - minus[[r]][-1L, , drop = FALSE]) / (2 * step)
      }
    }
    for (r in seq_len(nrow(cases))) for (h in 2:horizon) {
      jh <- matrix(jac[[r]][h - 1L, , ], k, n_par)
      se[[r]][h, ] <- sqrt(pmax(rowSums((jh %*% v$vcov) * jh), 0))
    }
  }
  unlist(lapply(se, function(m) as.vector(t(m))), use.names = FALSE)
}

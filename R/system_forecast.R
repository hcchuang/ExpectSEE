.system_stage <- function(fit) {
  fit <- active_stage(fit)
  if (!inherits(fit, c("expectsee_sals", "expectsee_sfgals",
                       "expectsee_hdsals", "expectsee_hdsfgals"))) {
    stop("fit must be a SALS/SFGALS or HDSALS/HDSFGALS estimate", call. = FALSE)
  }
  fit
}

.system_coefficients <- function(fit, equations, terms, residuals = FALSE) {
  rows <- fit$coef
  required <- c("equation", "term", "estimate")
  if (!is.data.frame(rows) || !all(required %in% names(rows)) ||
      !is.numeric(rows$estimate) || any(!is.finite(rows$estimate)) ||
      anyNA(rows[, c("equation", "term")]) ||
      !setequal(rows$equation, equations) ||
      anyDuplicated(rows[, c("equation", "term")])) {
    stop("invalid system coefficient table", call. = FALSE)
  }
  out <- matrix(0, length(equations), length(terms), dimnames = list(equations, terms))
  d <- matrix(0, length(equations), length(equations), dimnames = list(equations, equations))
  for (r in seq_len(nrow(rows))) {
    eq <- as.character(rows$equation[r])
    term <- as.character(rows$term[r])
    if (term == "(Intercept)") term <- "const"
    if (term %in% terms) {
      out[eq, term] <- rows$estimate[r]
    } else if (residuals && term %in% paste0("resid_", equations)) {
      source <- sub("^resid_", "", term)
      if (source == eq) stop("own-equation residual is not permitted", call. = FALSE)
      d[eq, source] <- rows$estimate[r]
    } else {
      stop("unsupported system term: ", term, call. = FALSE)
    }
  }
  if (any(vapply(equations, function(eq) {
    sum(rows$equation == eq & rows$term %in% c("const", "(Intercept)")) != 1L
  }, logical(1)))) stop("each equation must have exactly one intercept", call. = FALSE)
  list(coefficients = out, D = d)
}

#' Reduced system map for SEE forecasts and scenario responses
#'
#' Reconstructs both lag maps on a common basis and solves (I-D)G = C2-D*C1.
#' Invertibility is checked without a generalized-inverse fallback.
see_system_map <- function(fit, rcond_tol = 1e-8) {
  fit <- .system_stage(fit)
  stage2 <- inherits(fit, c("expectsee_sfgals", "expectsee_hdsfgals"))
  first <- if (stage2) fit$first_stage else fit
  if (!inherits(first, c("expectsee_sals", "expectsee_hdsals"))) {
    stop("fit has no valid first stage", call. = FALSE)
  }
  equations <- colnames(first$y)
  if (!length(equations) || anyNA(equations) || any(!nzchar(equations)) ||
      anyDuplicated(equations) || !setequal(colnames(fit$y), equations)) {
    stop("system equation names must be present, unique and aligned", call. = FALSE)
  }
  q <- first$q
  if (length(q) != 1L || !is.numeric(q) || !is.finite(q) || q < 1 || q != as.integer(q)) {
    stop("fit$q must be a positive integer", call. = FALSE)
  }
  if (!identical(fit$q, first$q) || !identical(fit$spec, first$spec) ||
      !identical(fit$tau, first$tau)) stop("stages have incompatible specifications", call. = FALSE)
  if (stage2 && (nrow(fit$y) != nrow(first$y) ||
      !isTRUE(all.equal(unname(as.matrix(fit$y)[, equations, drop = FALSE]),
                       unname(as.matrix(first$y)[, equations, drop = FALSE]), tolerance = 0)) ||
      !identical(fit$dates, first$dates))) {
    stop("first and second stage samples must be aligned", call. = FALSE)
  }
  check_tau(fit$tau)
  if (length(rcond_tol) != 1L || !is.numeric(rcond_tol) || !is.finite(rcond_tol) ||
      rcond_tol <= 0 || rcond_tol >= 1) stop("rcond_tol must be in (0, 1)", call. = FALSE)
  dummy <- matrix(0, q + 1L, length(equations), dimnames = list(NULL, equations))
  basis <- see_build_design(dummy, spec = first$spec, q = q)$colmap
  c1 <- .system_coefficients(first, equations, basis$term)$coefficients
  second <- .system_coefficients(fit, equations, basis$term, residuals = stage2)
  d <- second$D
  m <- diag(length(equations)) - d
  rc <- rcond(m)
  if (!is.finite(rc) || rc <= rcond_tol) {
    stop("I-D is singular or numerically ill-conditioned (rcond = ",
         format(rc), ")", call. = FALSE)
  }
  a <- solve(m)
  dimnames(a) <- dimnames(d)
  reduced <- solve(m, second$coefficients - d %*% c1)
  dimnames(reduced) <- dimnames(c1)
  if (any(!is.finite(c(a, reduced)))) stop("reduced system is not finite", call. = FALSE)
  reduced_residuals <- NULL
  if (!is.null(fit$resid)) {
    eta <- as.matrix(fit$resid)
    if (!is.numeric(eta) || !setequal(colnames(eta), equations) ||
        anyDuplicated(colnames(eta)) || nrow(eta) != nrow(fit$y) ||
        any(!is.finite(eta))) stop("invalid training residuals", call. = FALSE)
    reduced_residuals <- eta[, equations, drop = FALSE] %*% t(a)
  }
  rho <- max(Mod(eigen(d, only.values = TRUE)$values))
  structure(list(
    equations = equations, tau = fit$tau, spec = fit$spec, q = as.integer(q),
    estimator = sub("^expectsee_", "", class(fit)[1L]),
    coefficient_fit = fit$coefficient_fit %||% NA_character_,
    basis = basis, stage1_coefficients = c1, stage2_coefficients = second$coefficients,
    D = d, A = a, reduced_coefficients = reduced,
    reduced_residuals = reduced_residuals, history = first$y,
    estimation_end = if (length(first$dates)) as.character(utils::tail(first$dates, 1L)) else NA_character_,
    diagnostics = list(rcond = rc, inverse_norm_inf = norm(a, "I"),
                       spectral_radius_D = rho, neumann_series_converges = rho < 1),
    inference_available = FALSE
  ), class = "expectsee_system_map")
}

.system_history <- function(system, history = NULL) {
  history <- as.matrix(history %||% system$history)
  if (!is.numeric(history) || nrow(history) < system$q ||
      ncol(history) != length(system$equations) || any(!is.finite(history))) {
    stop("history must be finite numeric data with at least q rows and K columns", call. = FALSE)
  }
  if (is.null(colnames(history))) colnames(history) <- system$equations
  if (anyDuplicated(colnames(history)) || !setequal(colnames(history), system$equations)) {
    stop("history columns must match the system equations", call. = FALSE)
  }
  utils::tail(history[, system$equations, drop = FALSE], system$q)
}

.system_features <- function(system, history, columns = seq_len(nrow(system$basis))) {
  vapply(columns, function(i) {
    b <- system$basis[i, ]
    if (b$transform == "intercept") return(1)
    transform_lag(history[nrow(history) - b$lag + 1L, b$source_equation], b$transform)
  }, numeric(1))
}

.system_next <- function(system, history, correction) {
  active <- which(colSums(abs(system$reduced_coefficients)) > 0)
  value <- if (length(active)) {
    as.vector(system$reduced_coefficients[, active, drop = FALSE] %*%
                .system_features(system, history, active)) + correction
  } else correction
  if (any(!is.finite(value))) stop("recursive system path is not finite", call. = FALSE)
  stats::setNames(value, system$equations)
}

.system_vector <- function(x, equations, label) {
  if (!is.numeric(x) || !is.null(dim(x)) || length(x) != length(equations) || any(!is.finite(x))) {
    stop(label, " must be one finite numeric value per equation", call. = FALSE)
  }
  if (!is.null(names(x))) {
    if (anyDuplicated(names(x)) || !setequal(names(x), equations)) {
      stop(label, " names must match the equations", call. = FALSE)
    }
    x <- x[equations]
  }
  stats::setNames(as.numeric(x), equations)
}

.system_calibration <- function(system, calibration, innovation_expectile) {
  if (calibration == "provided") {
    return(.system_vector(innovation_expectile, system$equations, "innovation_expectile"))
  }
  if (!is.null(innovation_expectile)) {
    stop("use calibration = 'provided' with innovation_expectile", call. = FALSE)
  }
  u <- system$reduced_residuals
  if (is.null(u) || !nrow(u)) stop("in-sample calibration needs training residuals", call. = FALSE)
  correction <- apply(u, 2L, function(x) {
    f <- fit_expectile(x, matrix(1, length(x), 1L), system$tau)
    if (!f$converged) stop("innovation expectile calibration did not converge", call. = FALSE)
    f$beta[1L]
  })
  .system_vector(correction, system$equations, "calibration")
}

#' One-step SEE forecast using the reduced system
#'
#' Requires the conditional expectile of A*eta to equal the specified constant
#' correction. The default estimates it only from the fitted training residuals.
#' The supplied fit and history must both be available at the forecast origin.
see_forecast <- function(fit, history = NULL, horizon = 1L,
                         calibration = c("in_sample", "provided"),
                         innovation_expectile = NULL, rcond_tol = 1e-8) {
  if (length(horizon) != 1L || !is.numeric(horizon) || !is.finite(horizon) || horizon != 1) {
    stop("only one-step conditional expectile forecasts are supported", call. = FALSE)
  }
  calibration <- match.arg(calibration)
  system <- see_system_map(fit, rcond_tol)
  state <- .system_history(system, history)
  correction <- .system_calibration(system, calibration, innovation_expectile)
  reduced <- .system_next(system, state, rep(0, length(system$equations)))
  forecast <- reduced + correction
  if (any(!is.finite(forecast))) stop("system forecast is not finite", call. = FALSE)
  origin <- if (!is.null(rownames(state))) utils::tail(rownames(state), 1L) else NA_character_
  out <- data.frame(
    equation = system$equations, forecast = unname(forecast),
    reduced_component = unname(reduced), innovation_expectile = unname(correction),
    tau = system$tau, horizon = 1L, estimator = system$estimator,
    coefficient_fit = system$coefficient_fit, calibration = calibration,
    forecast_origin = origin, estimation_end = system$estimation_end,
    rcond = system$diagnostics$rcond, inverse_norm_inf = system$diagnostics$inverse_norm_inf,
    forecast_type = "one_step_constant_innovation_expectile", inference_available = FALSE,
    stringsAsFactors = FALSE
  )
  attr(out, "assumption") <- "conditional expectile of A*eta equals the correction at the forecast origin"
  class(out) <- c("expectsee_forecast", "data.frame")
  out
}


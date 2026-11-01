coef_vector_for_equation <- function(fit, equation, terms) {
  rows <- fit$coef[fit$coef$equation == equation, , drop = FALSE]
  if (!nrow(rows)) stop("coefficient table has no row for equation: ", equation, call. = FALSE)
  beta <- stats::setNames(rows$estimate, rows$term)
  missing <- setdiff(terms, names(beta))
  if (length(missing)) {
    stop("coefficient table is missing terms for equation ", equation, ": ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  beta[terms]
}

prediction_dates <- function(x) {
  out <- tryCatch(as.Date(x), error = function(e) NULL)
  if (!is.null(out) && !anyNA(out)) out else x
}

predict_stage1_object <- function(fit, y) {
  d <- see_build_design(y, spec = fit$spec, q = fit$q, include_intercept = TRUE)
  equations <- colnames(d$y)
  fitted <- resid <- matrix(NA_real_, nrow(d$y), ncol(d$y), dimnames = dimnames(d$y))

  for (i in seq_along(equations)) {
    terms <- fit$term_names[[i]]
    x <- matrix(NA_real_, nrow(d$y), length(terms), dimnames = list(NULL, terms))
    for (tm in terms) {
      if (tm %in% colnames(d$x)) {
        x[, tm] <- d$x[, tm]
      } else if (tm %in% c("const", "(Intercept)")) {
        x[, tm] <- 1
      } else {
        stop("prediction design is missing term: ", tm, call. = FALSE)
      }
    }
    beta <- coef_vector_for_equation(fit, equations[i], terms)
    fitted[, i] <- as.numeric(x %*% beta)
    resid[, i] <- d$y[, i] - fitted[, i]
  }

  structure(
    list(
      y = d$y, dates = prediction_dates(d$dates), fitted = fitted, resid = resid,
      design = d, spec = fit$spec, q = fit$q, tau = fit$tau,
      stage = if (inherits(fit, "expectsee_hdsals")) "hdsals" else "sals"
    ),
    class = "expectsee_prediction"
  )
}

predict_stage2_object <- function(fit, y) {
  first <- predict_stage1_object(fit$first_stage, y)
  equations <- colnames(first$y)
  fitted <- resid <- matrix(NA_real_, nrow(first$y), ncol(first$y), dimnames = dimnames(first$y))

  for (i in seq_along(equations)) {
    terms <- fit$term_names[[i]]
    x <- matrix(NA_real_, nrow(first$y), length(terms), dimnames = list(NULL, terms))
    for (tm in terms) {
      if (tm %in% colnames(first$design$x)) {
        x[, tm] <- first$design$x[, tm]
      } else if (tm %in% c("const", "(Intercept)")) {
        x[, tm] <- 1
      } else if (startsWith(tm, "resid_")) {
        source <- sub("^resid_", "", tm)
        if (!source %in% colnames(first$resid)) {
          stop("prediction residual matrix is missing source equation: ", source, call. = FALSE)
        }
        x[, tm] <- first$resid[, source]
      } else {
        stop("prediction design is missing term: ", tm, call. = FALSE)
      }
    }
    beta <- coef_vector_for_equation(fit, equations[i], terms)
    fitted[, i] <- as.numeric(x %*% beta)
    resid[, i] <- first$y[, i] - fitted[, i]
  }

  structure(
    list(
      y = first$y, dates = first$dates, fitted = fitted, resid = resid,
      first_stage = first, spec = fit$spec, q = fit$q, tau = fit$tau,
      stage = if (inherits(fit, "expectsee_hdsfgals")) "hdsfgals" else "sfgals"
    ),
    class = "expectsee_prediction"
  )
}

see_predict <- function(fit, y, stage = NULL) {
  obj <- active_stage(fit, stage)
  if (inherits(obj, c("expectsee_sals", "expectsee_hdsals"))) {
    return(predict_stage1_object(obj, y))
  }
  if (inherits(obj, c("expectsee_sfgals", "expectsee_hdsfgals"))) {
    return(predict_stage2_object(obj, y))
  }
  stop("fit must be an ExpectSEE stage fit or expectsee wrapper", call. = FALSE)
}

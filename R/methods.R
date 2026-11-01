print.expectsee <- function(x, ...) {
  cat("ExpectSEE fit\n")
  cat("  estimator:", x$estimator, "\n")
  if (!is.null(x$coefficient_fit) && !is.na(x$coefficient_fit)) {
    cat("  coefficient fit:", x$coefficient_fit_label, "\n")
  }
  cat("  spec:", x$spec, " tau:", x$tau, " q:", x$q, "\n")
  cat("  equations:", ncol(x$fitted), " observations:", nrow(x$fitted), "\n")
  invisible(x)
}

summary.expectsee <- function(object, ...) {
  list(
    call = object$call,
    coef = object$coef,
    gof = object$gof,
    diagnostics = object$diagnostics
  )
}

see_diagnostics <- function(fit, stage = NULL) {
  active_stage(fit, stage)$diagnostics
}

coef.expectsee <- function(object, stage = NULL, ...) {
  if (is.null(stage)) return(object$coef)
  object[[stage]]$coef
}

fitted.expectsee <- function(object, stage = NULL, ...) {
  if (is.null(stage)) return(object$fitted)
  object[[stage]]$fitted
}

residuals.expectsee <- function(object, stage = NULL, ...) {
  if (is.null(stage)) return(object$resid)
  object[[stage]]$resid
}

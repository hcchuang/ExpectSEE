active_stage <- function(fit, stage = NULL) {
  if (inherits(fit, "expectsee")) {
    if (is.null(stage)) {
      stage <- if (!is.null(fit$hdsfgals)) "hdsfgals" else if (!is.null(fit$sfgals)) "sfgals" else if (!is.null(fit$hdsals)) "hdsals" else "sals"
    }
    return(fit[[stage]])
  }
  fit
}

fixed_inference_object <- function(fit, stage = NULL) {
  obj <- active_stage(fit, stage)
  if (inherits(obj, c("expectsee_hdsals", "expectsee_hdsfgals"))) {
    stop("Wald inference is unavailable for high-dimensional Lasso fits.", call. = FALSE)
  }
  if (!inherits(obj, c("expectsee_sals", "expectsee_sfgals"))) {
    stop("fit must be a fixed-dimensional SALS or SFGALS object", call. = FALSE)
  }
  obj
}

equation_index <- function(obj, equation) {
  eqs <- colnames(obj$y)
  i <- if (is.character(equation)) match(equation, eqs) else as.integer(equation)
  if (length(i) != 1 || is.na(i) || i < 1 || i > length(eqs)) {
    stop("equation not found", call. = FALSE)
  }
  i
}

align_restrictions <- function(rmat, terms) {
  if (is.null(dim(rmat))) {
    names_r <- names(rmat)
    rmat <- matrix(rmat, nrow = 1)
    if (!is.null(names_r)) colnames(rmat) <- names_r
  }
  rmat <- as.matrix(rmat)
  if (!is.null(colnames(rmat))) {
    missing <- setdiff(colnames(rmat), terms)
    if (length(missing)) stop("restriction terms not found: ", paste(missing, collapse = ", "))
    out <- matrix(0, nrow(rmat), length(terms), dimnames = list(rownames(rmat), terms))
    out[, colnames(rmat)] <- rmat
    return(out)
  }
  if (ncol(rmat) != length(terms)) stop("R must have one column per coefficient")
  colnames(rmat) <- terms
  rmat
}

see_linear_test <- function(fit, equation, R, rhs = 0, stage = NULL,
                            se = c("generated", "sandwich")) {
  se <- match.arg(se)
  obj <- fixed_inference_object(fit, stage)
  i <- equation_index(obj, equation)
  terms <- obj$term_names[[i]]
  rmat <- align_restrictions(R, terms)
  rhs <- rep(rhs, length.out = nrow(rmat))
  if (length(rhs) != nrow(rmat) || any(!is.finite(rhs))) stop("invalid rhs")

  beta <- stats::setNames(obj$fits[[i]]$beta, terms)
  use_generated <- se == "generated" && !is.null(obj$vcov_generated)
  vc <- if (use_generated) obj$vcov_generated[[i]] else obj$vcov_sandwich[[i]]
  middle <- rmat %*% vc %*% t(rmat)
  if (any(!is.finite(middle))) stop("restriction covariance has invalid values")
  df <- qr(middle)$rank
  if (!df) stop("restriction covariance has zero rank")
  diff <- as.numeric(rmat %*% beta - rhs)
  statistic <- as.numeric(t(diff) %*% safe_inv(middle) %*% diff)

  data.frame(
    equation = colnames(obj$y)[i],
    se = if (use_generated) "generated" else "sandwich",
    statistic = statistic, df = df,
    p = stats::pchisq(statistic, df, lower.tail = FALSE),
    row.names = NULL
  )
}

see_wald_test <- function(fit, equation, terms = NULL, pattern = NULL,
                          stage = NULL, se = c("generated", "sandwich"), rhs = 0) {
  se <- match.arg(se)
  obj <- fixed_inference_object(fit, stage)
  i <- equation_index(obj, equation)
  names_i <- obj$term_names[[i]]
  idx <- if (!is.null(terms)) match(terms, names_i) else grep(pattern %||% ".", names_i)
  idx <- idx[!is.na(idx)]
  if (!length(idx)) stop("no tested terms found")
  rmat <- diag(length(names_i))[idx, , drop = FALSE]
  out <- see_linear_test(obj, equation, rmat, rhs, se = se)
  out$terms <- paste(names_i[idx], collapse = ", ")
  out
}

spec_restrictions <- function(obj, i, against) {
  map <- if (inherits(obj, "expectsee_sals")) {
    obj$term_maps[[i]]
  } else {
    obj$first_stage$term_maps[[i]]
  }
  terms <- obj$term_names[[i]]
  rows <- list()

  if (against %in% c("ABN", "ABP", "AB")) {
    transforms <- switch(
      against,
      ABN = "plus2",
      ABP = "minus2",
      AB = c("plus2", "minus2")
    )
    tested <- map$term[map$transform %in% transforms]
    rows <- lapply(tested, function(term) stats::setNames(1, term))
  } else {
    groups <- unique(map[map$transform %in% c("plus", "minus"),
                         c("source_equation", "source_index", "lag"), drop = FALSE])
    rows <- lapply(seq_len(nrow(groups)), function(g) {
      take <- !is.na(map$source_index) & !is.na(map$lag) &
        map$source_index == groups$source_index[g] & map$lag == groups$lag[g] &
        map$transform %in% c("plus", "minus")
      pair <- map$term[take]
      if (length(pair) != 2) return(NULL)
      stats::setNames(rep(1, 2), pair)
    })
    rows <- Filter(Negate(is.null), rows)
  }
  if (!length(rows)) stop("no specification restrictions found")
  rmat <- matrix(0, length(rows), length(terms), dimnames = list(NULL, terms))
  for (r in seq_along(rows)) rmat[r, names(rows[[r]])] <- rows[[r]]
  rmat
}

see_spec_test <- function(fit, against = c("ABN", "ABP", "AB", "SQ"), equation = NULL,
                          stage = NULL, se = c("generated", "sandwich")) {
  against <- match.arg(against)
  se <- match.arg(se)
  obj <- fixed_inference_object(fit, stage)
  if (!identical(obj$spec, "ABSQ")) stop("specification tests require an ABSQ fit")
  eqs <- colnames(obj$y)
  equations <- equation %||% eqs
  out <- lapply(equations, function(eq) {
    i <- equation_index(obj, eq)
    rmat <- spec_restrictions(obj, i, against)
    ans <- see_linear_test(obj, eq, rmat, se = se)
    ans$comparison <- paste("ABSQ vs", against)
    ans$terms <- paste(colnames(rmat)[colSums(abs(rmat)) > 0], collapse = ", ")
    ans
  })
  do.call(rbind, out)
}

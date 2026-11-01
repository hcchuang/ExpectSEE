see_specs <- function() c("AB", "SQ", "ABP", "ABN", "ABSQ")

normalize_spec <- function(spec) {
  spec <- toupper(spec)
  aliases <- c(ABS = "AB", PSQ = "ABP", NSQ = "ABN", BOTH = "ABSQ")
  if (spec %in% names(aliases)) spec <- aliases[[spec]]
  if (!spec %in% see_specs()) stop("spec must be one of: ", paste(see_specs(), collapse = ", "))
  spec
}

`%||%` <- function(x, y) if (is.null(x)) y else x

spec_transforms <- function(spec) {
  switch(
    normalize_spec(spec),
    AB = c("plus", "minus"),
    SQ = c("raw", "plus2", "minus2"),
    ABP = c("plus", "minus", "plus2"),
    ABN = c("plus", "minus", "minus2"),
    ABSQ = c("plus", "plus2", "minus", "minus2")
  )
}

transform_lag <- function(a, transform) {
  switch(
    transform,
    raw = a,
    plus = pmax(a, 0),
    minus = pmax(-a, 0),
    plus2 = pmax(a, 0)^2,
    minus2 = pmax(-a, 0)^2
  )
}

see_build_design <- function(y, spec = "ABSQ", q = 1, equation = NULL,
                             own_only = FALSE, include_intercept = TRUE) {
  spec <- normalize_spec(spec)
  y <- as.data.frame(y)
  if (is.null(colnames(y))) colnames(y) <- paste0("eq", seq_len(ncol(y)))
  if (q < 1 || q != as.integer(q)) stop("q must be a positive integer")
  if (nrow(y) <= q) stop("not enough rows for q lags")

  rows <- (q + 1):nrow(y)
  x_parts <- list()
  meta <- list()

  if (include_intercept) {
    x_parts[[length(x_parts) + 1]] <- matrix(rep(1, length(rows)), ncol = 1,
                                             dimnames = list(NULL, "const"))
    meta[[length(meta) + 1]] <- data.frame(
      term = "const", source_equation = NA_character_, source_index = NA_integer_,
      lag = NA_integer_, transform = "intercept", stringsAsFactors = FALSE
    )
  }

  for (lag in seq_len(q)) {
    lag_y <- as.matrix(y[rows - lag, , drop = FALSE])
    for (tr in spec_transforms(spec)) {
      val <- transform_lag(lag_y, tr)
      terms <- paste0("L", lag, "_", colnames(y), "_", tr)
      colnames(val) <- terms
      x_parts[[length(x_parts) + 1]] <- val
      meta[[length(meta) + 1]] <- data.frame(
        term = terms,
        source_equation = colnames(y),
        source_index = seq_along(y),
        lag = lag,
        transform = tr,
        stringsAsFactors = FALSE
      )
    }
  }

  x <- do.call(cbind, x_parts)
  colmap <- do.call(rbind, meta)
  rownames(colmap) <- NULL

  if (own_only) {
    if (is.null(equation)) stop("equation is required when own_only = TRUE")
    eq_id <- if (is.character(equation)) match(equation, colnames(y)) else equation
    keep <- is.na(colmap$source_index) | colmap$source_index == eq_id
    x <- x[, keep, drop = FALSE]
    colmap <- colmap[keep, , drop = FALSE]
  }

  list(
    x = x,
    y = as.matrix(y[rows, , drop = FALSE]),
    dates = rownames(y)[rows],
    colmap = colmap,
    spec = spec,
    q = q
  )
}

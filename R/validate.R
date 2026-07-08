# Author: Alberto Quaini

.as_numeric_matrix <- function(x, name = "x") {
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  } else {
    x <- as.matrix(x)
  }
  storage.mode(x) <- "double"
  if (!is.numeric(x) || length(dim(x)) != 2L) {
    stop(name, " must be coercible to a numeric matrix.", call. = FALSE)
  }
  if (nrow(x) < 1L || ncol(x) < 2L) {
    stop(name, " must have at least one row and two columns.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(name, " must contain only finite, non-missing values.", call. = FALSE)
  }
  x
}

.check_scalar_integer <- function(x, name, lower = -Inf, upper = Inf) {
  if (length(x) != 1L || is.na(x) || x != as.integer(x) || x < lower || x > upper) {
    stop(name, " must be an integer between ", lower, " and ", upper, ".", call. = FALSE)
  }
  as.integer(x)
}

.unit_vector <- function(x, name = "x") {
  x <- as.numeric(x)
  if (anyNA(x) || any(!is.finite(x))) {
    stop(name, " must contain only finite, non-missing values.", call. = FALSE)
  }
  nrm <- sqrt(sum(x^2))
  if (!is.finite(nrm) || nrm <= sqrt(.Machine$double.eps)) {
    stop(name, " must have positive Euclidean norm.", call. = FALSE)
  }
  x / nrm
}

.validate_rank <- function(p, max_rank) {
  .check_scalar_integer(p, "p", lower = 0L, upper = max_rank)
}

.check_fit <- function(x) {
  if (!inherits(x, "agca_fit")) {
    stop("x must be an object returned by agca() or agca_fit_directions().", call. = FALSE)
  }
  invisible(x)
}

#' Rank-Pareto marginal standardization
#'
#' Transforms each margin to empirical standard Pareto scores using
#' \eqn{(n + 1) / (n + 1 - rank)}. Larger observations are treated as more
#' extreme.
#'
#' @param x A numeric matrix or data frame.
#' @param ties_method Tie method passed to [base::rank()].
#' @return A numeric matrix with standard Pareto-like margins.
#' @export
rank_pareto <- function(x, ties_method = "average") {
  x <- .as_numeric_matrix(x)
  n <- nrow(x)
  out <- matrix(NA_real_, nrow = n, ncol = ncol(x), dimnames = dimnames(x))
  for (j in seq_len(ncol(x))) {
    ranks <- rank(x[, j], ties.method = ties_method)
    out[, j] <- (n + 1) / (n + 1 - ranks)
  }
  out
}

#' Pareto marginal standardization from supplied CDFs
#'
#' @param x A numeric matrix or data frame.
#' @param cdf A function applied to every margin, or a list of one CDF function
#'   per margin. Each CDF must return values in `[0, 1]`.
#' @param eps Tail clipping constant used to avoid zero and infinite values.
#' @return A numeric matrix with standard Pareto margins.
#' @export
pareto_from_cdf <- function(x, cdf, eps = 1e-12) {
  x <- .as_numeric_matrix(x)
  if (length(eps) != 1L || !is.finite(eps) || eps <= 0 || eps >= 0.5) {
    stop("eps must be a finite number in (0, 0.5).", call. = FALSE)
  }

  if (is.function(cdf)) {
    cdf <- rep(list(cdf), ncol(x))
  }
  if (!is.list(cdf) || length(cdf) != ncol(x) || !all(vapply(cdf, is.function, logical(1L)))) {
    stop("cdf must be a function or a list of one function per margin.", call. = FALSE)
  }

  out <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  for (j in seq_len(ncol(x))) {
    u <- cdf[[j]](x[, j])
    if (length(u) != nrow(x) || anyNA(u) || any(!is.finite(u)) || any(u < 0 | u > 1)) {
      stop("Each CDF must return finite values in [0, 1].", call. = FALSE)
    }
    u <- pmin(1 - eps, pmax(eps, u))
    out[, j] <- 1 / (1 - u)
  }
  out
}

#' Marginal standardization for AGCA
#'
#' @param x A numeric matrix or data frame.
#' @param margin Standardization method: `"rank_pareto"` (default),
#'   `"pareto"` for supplied CDFs, or `"none"` for already standardized data.
#' @param cdf Optional CDF function or list of CDF functions for
#'   `margin = "pareto"`.
#' @param ties_method Tie method used by [rank_pareto()].
#' @return A numeric matrix.
#' @export
agca_standardize <- function(x, margin = c("rank_pareto", "pareto", "none"),
                             cdf = NULL, ties_method = "average") {
  margin <- match.arg(margin)
  switch(
    margin,
    rank_pareto = rank_pareto(x, ties_method = ties_method),
    pareto = {
      if (is.null(cdf)) {
        stop("cdf must be supplied when margin = 'pareto'.", call. = FALSE)
      }
      pareto_from_cdf(x, cdf = cdf)
    },
    none = .as_numeric_matrix(x)
  )
}

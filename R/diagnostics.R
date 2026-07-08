#' Threshold stability diagnostics
#'
#' @param x Numeric matrix or data frame.
#' @param k Integer vector of top-k values.
#' @param ... Additional arguments passed to [agca()].
#' @return A data frame of rank summaries across thresholds.
#' @export
threshold_stability <- function(x, k, ...) {
  k <- as.integer(k)
  if (anyNA(k) || any(k < 1L)) {
    stop("k must contain positive integers.", call. = FALSE)
  }
  out <- lapply(k, function(kk) {
    fit <- agca(x, k = kk, bootstrap = NULL, ...)
    ans <- agca_rank_summary(fit)
    ans$k <- kk
    ans
  })
  do.call(rbind, out)
}

#' Anchor sensitivity diagnostics
#'
#' @param x Numeric matrix or data frame.
#' @param k Number of top radial observations.
#' @param anchors Character vector of anchors to compare.
#' @param ... Additional arguments passed to [agca()].
#' @return A data frame of rank summaries across anchors.
#' @export
anchor_sensitivity <- function(x, k, anchors = c("canonical", "frechet", "principal"), ...) {
  out <- lapply(anchors, function(a) {
    fit <- agca(x, k = k, anchor = a, bootstrap = NULL, ...)
    ans <- agca_rank_summary(fit)
    ans$anchor = a
    ans
  })
  do.call(rbind, out)
}

#' Functional approximation errors from angular reconstruction
#'
#' Computes mean angular functional values on the fitted and reconstructed
#' directions for a collection of portfolio weights.
#'
#' @param fit An object returned by [agca()] or [agca_fit_directions()].
#' @param weights A numeric vector or matrix. Rows are portfolios.
#' @param ranks Integer vector of AGCA ranks.
#' @param power Power applied to positive portfolio exposures.
#' @param cap Optional finite cap applied to the powered exposure.
#' @return A data frame with original, reconstructed, and relative errors.
#' @export
angular_functional_error <- function(fit, weights, ranks = fit$p,
                                     power = 1, cap = Inf) {
  .check_fit(fit)
  weights <- as.matrix(weights)
  storage.mode(weights) <- "double"
  if (ncol(weights) != fit$dimension) {
    if (nrow(weights) == fit$dimension) {
      weights <- t(weights)
    } else {
      stop("weights must have one column per AGCA dimension.", call. = FALSE)
    }
  }
  ranks <- as.integer(ranks)
  if (anyNA(ranks) || any(ranks < 0L | ranks > length(fit$eigenvalues))) {
    stop("ranks must contain integers between 0 and the maximum AGCA rank.", call. = FALSE)
  }
  if (length(power) != 1L || !is.finite(power) || power <= 0) {
    stop("power must be a positive finite number.", call. = FALSE)
  }
  if (length(cap) != 1L || is.na(cap) || cap <= 0) {
    stop("cap must be positive or Inf.", call. = FALSE)
  }

  portfolio_names <- rownames(weights)
  if (is.null(portfolio_names)) {
    portfolio_names <- paste0("portfolio_", seq_len(nrow(weights)))
  }
  core <- agca_functional_error_cpp(
    fit$g,
    fit$anchor_coordinate,
    fit$mu,
    fit$scores,
    fit$loadings,
    weights,
    ranks,
    power,
    cap
  )

  data.frame(
    rank = rep(ranks, each = nrow(weights)),
    portfolio = rep(portfolio_names, times = length(ranks)),
    original = rep(as.numeric(core$original), times = length(ranks)),
    reconstructed = as.vector(t(core$reconstructed)),
    relative_error = as.vector(t(core$relative_error)),
    row.names = NULL
  )
}

# Author: Alberto Quaini

#' Nonparametric AGCA bootstrap
#'
#' Resamples fitted angular directions and recomputes AGCA diagnostics.
#'
#' @param fit An object returned by [agca()] or [agca_fit_directions()].
#' @param B Number of bootstrap resamples.
#' @param ranks Integer ranks to summarize.
#' @param fixed_anchor Logical. If `TRUE`, keep the fitted anchor fixed.
#'   Otherwise refit the requested anchor type.
#' @param anchor Anchor used when `fixed_anchor = FALSE`.
#' @param seed Optional random seed.
#' @return An object of class `"agca_bootstrap"`.
#' @export
bootstrap_agca <- function(fit, B = 199L, ranks = NULL, fixed_anchor = TRUE,
                           anchor = "canonical", seed = NULL) {
  .check_fit(fit)
  B <- .check_scalar_integer(B, "B", lower = 1L)
  max_rank <- length(fit$eigenvalues)
  if (is.null(ranks)) {
    ranks <- 0:max_rank
  }
  ranks <- as.integer(ranks)
  if (anyNA(ranks) || any(ranks < 0L | ranks > max_rank)) {
    stop("ranks must contain integers between 0 and the maximum AGCA rank.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (isTRUE(fixed_anchor)) {
    anchor_type <- "fixed"
    anchor_vector <- numeric(0L)
  } else if (is.character(anchor)) {
    anchor_type <- match.arg(anchor, c("canonical", "frechet", "principal"))
    anchor_vector <- numeric(0L)
  } else if (is.numeric(anchor)) {
    anchor_type <- "numeric"
    anchor_vector <- .unit_vector(anchor, "anchor")
  } else {
    stop("anchor must be 'canonical', 'frechet', 'principal', or a numeric vector.", call. = FALSE)
  }

  out <- list(
    replicates = agca_bootstrap_cpp(
      fit$g,
      fit$mu,
      p = fit$p,
      B = B,
      ranks = ranks,
      fixed_anchor = isTRUE(fixed_anchor),
      anchor_type = anchor_type,
      anchor_vector = anchor_vector
    ),
    B = B,
    ranks = ranks,
    fixed_anchor = isTRUE(fixed_anchor)
  )
  class(out) <- "agca_bootstrap"
  out
}

#' Summarize AGCA bootstrap output
#'
#' @param object An object returned by [bootstrap_agca()].
#' @param probs Quantile probabilities.
#' @param ... Unused.
#' @return A data frame of bootstrap summaries by rank.
#' @export
summary.agca_bootstrap <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  if (!inherits(object, "agca_bootstrap")) {
    stop("object must be returned by bootstrap_agca().", call. = FALSE)
  }
  probs <- as.numeric(probs)
  dat <- object$replicates
  by_rank <- split(dat, dat$rank)
  out <- lapply(by_rank, function(z) {
    risk_q <- stats::quantile(z$residual_risk, probs = probs, names = FALSE)
    ave_q <- stats::quantile(z$variation_explained, probs = probs, names = FALSE)
    data.frame(
      rank = z$rank[1L],
      statistic = c(rep("residual_risk", length(probs)), rep("variation_explained", length(probs))),
      prob = rep(probs, 2L),
      value = c(risk_q, ave_q),
      row.names = NULL
    )
  })
  do.call(rbind, out)
}

#' @export
print.agca_bootstrap <- function(x, ...) {
  cat("AGCA bootstrap\n")
  cat("  resamples:", x$B, "\n")
  cat("  ranks:", paste(x$ranks, collapse = ", "), "\n")
  cat("  fixed anchor:", x$fixed_anchor, "\n")
  invisible(x)
}

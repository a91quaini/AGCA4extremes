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

  n <- nrow(fit$g)
  records <- vector("list", B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    boot_anchor <- if (isTRUE(fixed_anchor)) fit$mu else anchor
    boot_fit <- agca_fit_directions(
      fit$g[idx, , drop = FALSE],
      anchor = boot_anchor,
      p = fit$p,
      normalize = FALSE
    )
    risk <- agca_residual_risk(boot_fit, max_rank = max_rank)
    ave <- c(0, agca_variation_explained(boot_fit))
    records[[b]] <- data.frame(
      replicate = b,
      rank = ranks,
      residual_risk = unname(risk[as.character(ranks)]),
      variation_explained = unname(ave[ranks + 1L]),
      row.names = NULL
    )
  }

  out <- list(
    replicates = do.call(rbind, records),
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

# Author: Alberto Quaini

#' Plot AGCA output
#'
#' @param x An object returned by [agca()] or [agca_fit_directions()].
#' @param type Plot type: eigenvalue scree plot, cumulative variation, first two
#'   scores, or loadings.
#' @param p Component index used for the loadings plot.
#' @param ... Additional graphical arguments.
#' @return Invisibly returns `x`.
#' @export
plot.agca_fit <- function(x, type = c("variation", "scree", "scores", "loadings"),
                          p = 1L, ...) {
  .check_fit(x)
  type <- match.arg(type)
  if (type == "variation") {
    ave <- agca_variation_explained(x)
    graphics::plot(
      seq_along(ave), ave,
      type = "b", xlab = "Rank", ylab = "Anchored variation explained",
      ylim = c(0, 1), ...
    )
  } else if (type == "scree") {
    graphics::plot(
      seq_along(x$eigenvalues), x$eigenvalues,
      type = "b", xlab = "Component", ylab = "Eigenvalue", ...
    )
  } else if (type == "scores") {
    if (ncol(x$scores) < 2L) {
      stop("At least two components are needed for a score plot.", call. = FALSE)
    }
    graphics::plot(
      x$scores[, 1L], x$scores[, 2L],
      xlab = "AGC1 score", ylab = "AGC2 score", ...
    )
  } else {
    p <- .check_scalar_integer(p, "p", lower = 1L, upper = ncol(x$loadings))
    graphics::barplot(
      x$loadings[, p],
      ylab = paste0("AGC", p, " loading"),
      las = 2,
      ...
    )
  }
  invisible(x)
}

#' Plot AGCA bootstrap summaries
#'
#' @param x An object returned by [bootstrap_agca()].
#' @param statistic Statistic to plot.
#' @param ... Additional graphical arguments.
#' @return Invisibly returns `x`.
#' @export
plot.agca_bootstrap <- function(x, statistic = c("variation_explained", "residual_risk"), ...) {
  if (!inherits(x, "agca_bootstrap")) {
    stop("x must be returned by bootstrap_agca().", call. = FALSE)
  }
  statistic <- match.arg(statistic)
  dat <- x$replicates
  graphics::boxplot(
    dat[[statistic]] ~ dat$rank,
    xlab = "Rank",
    ylab = statistic,
    ...
  )
  invisible(x)
}

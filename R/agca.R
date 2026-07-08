#' Fit anchored geodesic component analysis
#'
#' `agca()` is the main user interface. It standardizes margins, extracts
#' large-radius observations, resolves the anchor, and fits anchored geodesic
#' components to the resulting angular directions.
#'
#' @param x Numeric matrix or data frame. Larger values are treated as more
#'   extreme in each margin.
#' @param k Number of largest radial observations to retain. Specify either
#'   `k` or `threshold`.
#' @param threshold Radial threshold for selecting extremes.
#' @param margin Marginal standardization method. The default `"rank_pareto"`
#'   uses empirical ranks. Use `"pareto"` with `cdf`, or `"none"` for already
#'   standardized observations.
#' @param cdf Optional CDF function or list of CDF functions for
#'   `margin = "pareto"`.
#' @param anchor `"canonical"`, `"frechet"`, `"principal"`, or a numeric anchor
#'   vector.
#' @param p Working reconstruction rank. Defaults to the full tangent rank.
#' @param decluster Optional. Use `NULL` for no declustering, `TRUE` for runs
#'   declustering with `run = 1`, a nonnegative integer run length, or a list
#'   with component `run`.
#' @param bootstrap Optional integer number of bootstrap resamples.
#' @param keep_data Logical. If `TRUE`, store the standardized data matrix in
#'   the returned object.
#' @param ties_method Tie method used by [rank_pareto()].
#' @param seed Optional random seed used when `bootstrap` is supplied.
#' @return An object of class `"agca_fit"`.
#' @export
#' @examples
#' data(ff_portfolio_losses)
#' x <- ff_portfolio_losses[seq_len(500), -1]
#' fit <- agca(x, k = 75, p = 3)
#' fit
#' agca_rank_summary(fit)
agca <- function(x, k = NULL, threshold = NULL,
                 margin = c("rank_pareto", "pareto", "none"),
                 cdf = NULL, anchor = "canonical", p = NULL,
                 decluster = NULL, bootstrap = NULL,
                 keep_data = FALSE, ties_method = "average", seed = NULL) {
  matched_call <- match.call()
  margin <- match.arg(margin)
  z <- agca_standardize(x, margin = margin, cdf = cdf, ties_method = ties_method)

  if (is.null(k) && is.null(threshold)) {
    stop("Specify k or threshold for the extreme radial sample.", call. = FALSE)
  }

  tail <- if (is.null(decluster)) {
    tail_directions(z, k = k, threshold = threshold)
  } else {
    run <- if (isTRUE(decluster)) {
      1L
    } else if (is.numeric(decluster)) {
      decluster
    } else if (is.list(decluster) && !is.null(decluster$run)) {
      decluster$run
    } else {
      stop("decluster must be NULL, TRUE, a run length, or a list with run.", call. = FALSE)
    }
    decluster_runs(z, k = k, threshold = threshold, run = run)
  }

  fit <- agca_fit_directions(tail$directions, anchor = anchor, p = p, normalize = FALSE)
  fit$call <- matched_call
  fit$margin <- margin
  fit$tail <- tail
  fit$declustered <- !is.null(decluster)
  fit$standardized_data <- if (isTRUE(keep_data)) z else NULL

  if (!is.null(bootstrap)) {
    bootstrap <- .check_scalar_integer(bootstrap, "bootstrap", lower = 1L)
    fit$bootstrap <- bootstrap_agca(fit, B = bootstrap, seed = seed)
  }

  fit
}

#' Fit AGCA to angular directions
#'
#' @param g Matrix of angular directions.
#' @param anchor `"canonical"`, `"frechet"`, `"principal"`, or a numeric anchor
#'   vector.
#' @param p Working reconstruction rank. Defaults to the full tangent rank.
#' @param normalize Logical. If `TRUE`, rows of `g` are normalized first.
#' @return An object of class `"agca_fit"`.
#' @export
agca_fit_directions <- function(g, anchor = "canonical", p = NULL, normalize = TRUE) {
  g <- .as_numeric_matrix(g, "g")
  mu <- resolve_anchor(anchor, if (isTRUE(normalize)) normalize_rows(g) else g)
  core <- agca_core_cpp(g, mu, normalize = isTRUE(normalize))
  max_rank <- length(core$eigenvalues)
  if (is.null(p)) {
    p <- max_rank
  }
  p <- .validate_rank(p, max_rank)

  names(core$eigenvalues) <- paste0("AGC", seq_along(core$eigenvalues))
  rownames(core$loadings) <- colnames(core$g)
  colnames(core$loadings) <- paste0("AGC", seq_len(ncol(core$loadings)))
  colnames(core$scores) <- colnames(core$loadings)

  out <- c(
    core,
    list(
      p = p,
      n_extreme = nrow(core$g),
      dimension = ncol(core$g),
      anchor = core$mu,
      call = match.call(),
      margin = "none",
      tail = NULL,
      declustered = FALSE,
      standardized_data = NULL,
      bootstrap = NULL
    )
  )
  class(out) <- "agca_fit"
  out
}

#' Reconstruct angular directions from leading AGCA components
#'
#' @param fit An object returned by [agca()] or [agca_fit_directions()].
#' @param p Reconstruction rank.
#' @return A matrix of reconstructed angular directions.
#' @export
agca_reconstruct <- function(fit, p = fit$p) {
  .check_fit(fit)
  p <- .validate_rank(p, length(fit$eigenvalues))
  out <- agca_reconstruct_cpp(
    fit$anchor_coordinate,
    fit$mu,
    fit$scores,
    fit$loadings,
    p
  )
  rownames(out) <- rownames(fit$g)
  colnames(out) <- colnames(fit$g)
  out
}

#' Residual risk by AGCA rank
#'
#' @param fit An object returned by [agca()] or [agca_fit_directions()].
#' @param max_rank Maximum rank to report.
#' @return A numeric vector indexed by ranks `0:max_rank`.
#' @export
agca_residual_risk <- function(fit, max_rank = length(fit$eigenvalues)) {
  .check_fit(fit)
  max_rank <- .validate_rank(max_rank, length(fit$eigenvalues))
  out <- agca_residual_risk_cpp(fit$u, fit$scores, fit$loadings, max_rank)
  names(out) <- as.character(0:max_rank)
  out
}

#' Anchored variation explained
#'
#' @param fit An object returned by [agca()] or [agca_fit_directions()].
#' @return Cumulative anchored variation explained by each rank.
#' @export
agca_variation_explained <- function(fit) {
  .check_fit(fit)
  total <- sum(fit$eigenvalues)
  if (!is.finite(total) || total <= 0) {
    return(rep(NA_real_, length(fit$eigenvalues)))
  }
  out <- cumsum(fit$eigenvalues) / total
  names(out) <- names(fit$eigenvalues)
  out
}

#' AGCA rank summary
#'
#' @param fit An object returned by [agca()] or [agca_fit_directions()].
#' @return A data frame with rank, residual risk, and variation explained.
#' @export
agca_rank_summary <- function(fit) {
  .check_fit(fit)
  max_rank <- length(fit$eigenvalues)
  data.frame(
    rank = 0:max_rank,
    residual_risk = agca_residual_risk(fit, max_rank = max_rank),
    variation_explained = c(0, agca_variation_explained(fit)),
    row.names = NULL
  )
}

#' @export
print.agca_fit <- function(x, ...) {
  .check_fit(x)
  cat("Anchored geodesic component analysis\n")
  cat("  extremes:", x$n_extreme, "\n")
  cat("  dimension:", x$dimension, "\n")
  cat("  rank:", x$p, "\n")
  cat("  margin:", x$margin, "\n")
  cat("  declustered:", x$declustered, "\n")
  ave <- agca_variation_explained(x)
  shown <- seq_len(min(5L, length(ave)))
  tab <- data.frame(
    component = names(x$eigenvalues)[shown],
    eigenvalue = unname(x$eigenvalues[shown]),
    variation_explained = unname(ave[shown])
  )
  print(tab, row.names = FALSE, digits = 4)
  invisible(x)
}

#' @export
summary.agca_fit <- function(object, ...) {
  .check_fit(object)
  out <- list(
    call = object$call,
    n_extreme = object$n_extreme,
    dimension = object$dimension,
    p = object$p,
    margin = object$margin,
    declustered = object$declustered,
    rank_summary = agca_rank_summary(object)
  )
  class(out) <- "summary_agca_fit"
  out
}

#' @export
print.summary_agca_fit <- function(x, ...) {
  cat("AGCA summary\n")
  cat("  extremes:", x$n_extreme, "\n")
  cat("  dimension:", x$dimension, "\n")
  cat("  rank:", x$p, "\n\n")
  print(x$rank_summary, row.names = FALSE, digits = 4)
  invisible(x)
}

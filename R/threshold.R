#' Top-k angular directions
#'
#' @param x A numeric matrix of standardized observations.
#' @param k Number of largest radii to retain. Specify either `k` or
#'   `threshold`.
#' @param threshold Radial threshold. Observations with radius greater than
#'   `threshold` are retained.
#' @return A list containing angular directions, radii, selected indices, and
#'   the threshold.
#' @export
tail_directions <- function(x, k = NULL, threshold = NULL) {
  x <- .as_numeric_matrix(x)

  if (!is.null(k) && !is.null(threshold)) {
    stop("Specify either k or threshold, not both.", call. = FALSE)
  }
  if (is.null(k) && is.null(threshold)) {
    stop("Specify k or threshold.", call. = FALSE)
  }

  if (!is.null(k)) {
    k <- .check_scalar_integer(k, "k", lower = 1L, upper = nrow(x))
    out <- agca_tail_directions_cpp(x, k = k, threshold = NA_real_, use_k = TRUE)
  } else {
    if (length(threshold) != 1L || !is.finite(threshold)) {
      stop("threshold must be a finite number.", call. = FALSE)
    }
    out <- agca_tail_directions_cpp(x, k = 0L, threshold = threshold, use_k = FALSE)
  }

  colnames(out$directions) <- colnames(x)
  out
}

#' Runs declustering for radial extremes
#'
#' Exceedances are split into clusters separated by more than `run` consecutive
#' non-exceedances. The representative of each cluster is the observation with
#' the largest radius.
#'
#' @param x A numeric matrix of standardized observations.
#' @param k,threshold Top-k count or radial threshold used to define
#'   exceedances.
#' @param run Nonnegative run length.
#' @return A list like [tail_directions()], with one index per cluster.
#' @export
decluster_runs <- function(x, k = NULL, threshold = NULL, run = 1L) {
  x <- .as_numeric_matrix(x)
  run <- .check_scalar_integer(run, "run", lower = 0L)

  if (!is.null(k) && !is.null(threshold)) {
    stop("Specify either k or threshold, not both.", call. = FALSE)
  }
  if (is.null(k) && is.null(threshold)) {
    stop("Specify k or threshold.", call. = FALSE)
  }

  if (!is.null(k)) {
    k <- .check_scalar_integer(k, "k", lower = 1L, upper = nrow(x))
    out <- agca_decluster_runs_cpp(x, k = k, threshold = NA_real_, use_k = TRUE, run = run)
  } else {
    if (length(threshold) != 1L || !is.finite(threshold)) {
      stop("threshold must be a finite number.", call. = FALSE)
    }
    out <- agca_decluster_runs_cpp(x, k = 0L, threshold = threshold, use_k = FALSE, run = run)
  }
  colnames(out$directions) <- colnames(x)
  out
}

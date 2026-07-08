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
  radius <- row_norms(x)

  if (!is.null(k) && !is.null(threshold)) {
    stop("Specify either k or threshold, not both.", call. = FALSE)
  }
  if (is.null(k) && is.null(threshold)) {
    stop("Specify k or threshold.", call. = FALSE)
  }

  if (!is.null(k)) {
    k <- .check_scalar_integer(k, "k", lower = 1L, upper = nrow(x))
    index <- order(radius, decreasing = TRUE)[seq_len(k)]
    threshold <- min(radius[index])
  } else {
    if (length(threshold) != 1L || !is.finite(threshold)) {
      stop("threshold must be a finite number.", call. = FALSE)
    }
    index <- which(radius > threshold)
  }

  if (length(index) == 0L) {
    stop("No observations exceed the threshold.", call. = FALSE)
  }

  list(
    directions = normalize_rows(x[index, , drop = FALSE]),
    radius = radius[index],
    index = index,
    threshold = threshold,
    all_radius = radius
  )
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
  radius <- row_norms(x)
  tail <- tail_directions(x, k = k, threshold = threshold)
  exceed <- rep(FALSE, nrow(x))
  exceed[tail$index] <- TRUE
  exceed_index <- which(exceed)

  clusters <- list()
  current <- exceed_index[1L]
  for (idx in exceed_index[-1L]) {
    gap <- idx - current[length(current)] - 1L
    if (gap <= run) {
      current <- c(current, idx)
    } else {
      clusters[[length(clusters) + 1L]] <- current
      current <- idx
    }
  }
  clusters[[length(clusters) + 1L]] <- current

  index <- vapply(
    clusters,
    function(idx) idx[which.max(radius[idx])],
    integer(1L)
  )

  list(
    directions = normalize_rows(x[index, , drop = FALSE]),
    radius = radius[index],
    index = index,
    threshold = tail$threshold,
    all_radius = radius,
    clusters = clusters,
    run = run
  )
}

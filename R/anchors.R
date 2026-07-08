#' Canonical anchor
#'
#' @param d Ambient dimension.
#' @return The balanced direction \eqn{d^{-1/2}(1,\ldots,1)}.
#' @export
canonical_anchor <- function(d) {
  .check_scalar_integer(d, "d", lower = 2L)
  rep(1 / sqrt(d), d)
}

#' Principal angular anchor
#'
#' @param g Matrix of angular directions.
#' @param normalize Logical. If `TRUE`, rows are normalized first.
#' @return The leading eigenvector of the angular second-moment matrix,
#'   oriented to have positive sum.
#' @export
principal_anchor <- function(g, normalize = TRUE) {
  g <- .as_numeric_matrix(g, "g")
  if (isTRUE(normalize)) {
    g <- normalize_rows(g)
  }
  eig <- eigen(crossprod(g) / nrow(g), symmetric = TRUE)
  mu <- eig$vectors[, 1L]
  if (sum(mu) < 0) {
    mu <- -mu
  }
  .unit_vector(mu, "principal anchor")
}

#' Spherical Frechet anchor
#'
#' Computes a Karcher-mean approximation to the spherical Frechet mean.
#'
#' @param g Matrix of angular directions.
#' @param normalize Logical. If `TRUE`, rows are normalized first.
#' @param max_iter Maximum number of iterations.
#' @param tol Convergence tolerance for the tangent update norm.
#' @return A unit vector.
#' @export
frechet_anchor <- function(g, normalize = TRUE, max_iter = 100L, tol = 1e-10) {
  g <- .as_numeric_matrix(g, "g")
  if (isTRUE(normalize)) {
    g <- normalize_rows(g)
  }
  max_iter <- .check_scalar_integer(max_iter, "max_iter", lower = 1L)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be a positive finite number.", call. = FALSE)
  }

  mu <- .unit_vector(colMeans(g), "initial Frechet anchor")
  for (iter in seq_len(max_iter)) {
    inner <- as.numeric(g %*% mu)
    inner <- pmin(1, pmax(-1, inner))
    theta <- acos(inner)
    scale <- rep(1, length(theta))
    nz <- theta > sqrt(.Machine$double.eps)
    scale[nz] <- theta[nz] / sin(theta[nz])
    logs <- (g - tcrossprod(inner, mu)) * scale
    update <- colMeans(logs)
    update_norm <- sqrt(sum(update^2))
    if (update_norm < tol) {
      break
    }
    mu <- cos(update_norm) * mu + sin(update_norm) * update / update_norm
    mu <- .unit_vector(mu, "Frechet anchor")
  }
  if (sum(mu) < 0) {
    mu <- -mu
  }
  mu
}

resolve_anchor <- function(anchor, g) {
  d <- ncol(g)
  if (is.character(anchor)) {
    anchor <- match.arg(anchor, c("canonical", "frechet", "principal"))
    return(switch(
      anchor,
      canonical = canonical_anchor(d),
      frechet = frechet_anchor(g),
      principal = principal_anchor(g)
    ))
  }
  if (is.numeric(anchor)) {
    if (length(anchor) != d) {
      stop("Numeric anchor must have length equal to ncol(g).", call. = FALSE)
    }
    return(.unit_vector(anchor, "anchor"))
  }
  stop("anchor must be 'canonical', 'frechet', 'principal', or a numeric vector.", call. = FALSE)
}

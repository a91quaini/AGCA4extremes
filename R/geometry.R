#' Row Euclidean norms
#'
#' @param x A numeric matrix.
#' @return A numeric vector containing one Euclidean norm per row.
#' @export
row_norms <- function(x) {
  agca_row_norms_cpp(.as_numeric_matrix(x))
}

#' Normalize matrix rows
#'
#' @param x A numeric matrix with nonzero rows.
#' @return A numeric matrix whose rows have Euclidean norm one.
#' @export
normalize_rows <- function(x) {
  agca_normalize_rows_cpp(.as_numeric_matrix(x))
}

#' Spherical geodesic distance
#'
#' Computes great-circle distances on the unit sphere. If `y` has one row and
#' `x` has several rows, the single direction in `y` is recycled.
#'
#' @param x,y Numeric matrices with the same number of columns.
#' @param normalize Logical. If `TRUE`, rows are normalized before distances
#'   are computed.
#' @return A numeric vector of geodesic distances in radians.
#' @export
sphere_distance <- function(x, y, normalize = TRUE) {
  as.numeric(agca_geodesic_distance_cpp(
    .as_numeric_matrix(x, "x"),
    .as_numeric_matrix(y, "y"),
    normalize = isTRUE(normalize)
  ))
}

#' Anchored departures
#'
#' Computes the anchor coordinate and tangent departure
#' \eqn{u_\mu(g) = (I-\mu\mu^\top)g}.
#'
#' @param g A numeric matrix of directions.
#' @param mu Anchor direction.
#' @param normalize Logical. If `TRUE`, rows of `g` are normalized first.
#' @return A list with normalized directions, anchor, anchor coordinates, and
#'   anchored departures.
#' @export
anchored_departures <- function(g, mu, normalize = TRUE) {
  agca_departures_cpp(
    .as_numeric_matrix(g, "g"),
    .unit_vector(mu, "mu"),
    normalize = isTRUE(normalize)
  )
}

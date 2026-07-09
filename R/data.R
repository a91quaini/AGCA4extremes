# Author: Alberto Quaini

#' Simulated 10-dimensional AGCA example
#'
#' A package-owned simulated heavy-tailed sample from the 10-dimensional design
#' used in *Anchored Geodesic Analysis for Multivariate Extremes*. Variables
#' `X1`--`X8` share a low-dimensional logistic-block extremal mechanism.
#' Variables `X9` and `X10` contain independent Pareto sources, creating
#' near-axis extreme regimes alongside the shared low-rank angular structure.
#'
#' @format A data frame with 10,000 rows and 11 columns. Columns `X1`, ...,
#' `X10` are positive heavy-tailed observations. Column `regime` is a latent
#' factor identifying the dominant source for the observation: shared low-rank,
#' axis 9, or axis 10.
#' @source Simulated by `data-raw/simulate_data.R` using
#' [simulate_agca_10d()].
"agca_10d_simulation"

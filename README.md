# AGCA4extremes

Author: Alberto Quaini

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/AGCA4extremes)](https://CRAN.R-project.org/package=AGCA4extremes)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/a91quaini/AGCA4extremes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/a91quaini/AGCA4extremes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`AGCA4extremes` implements anchored geodesic component analysis (AGCA) for
multivariate extremes. AGCA summarizes how extreme angular profiles vary around
a benchmark direction, most often the canonical balanced-dependence anchor.

The package provides:

- rank-Pareto marginal standardization;
- Pareto standardization from user-supplied CDFs;
- top-`k` angular extraction and optional runs declustering;
- canonical, spherical Frechet, and principal anchors;
- AGCA scores, loadings, eigenvalues, explained variation, reconstruction, and
  residual risk;
- bootstrap, threshold, and anchor-sensitivity diagnostics;
- base R plotting methods and functional approximation summaries.

The numerical core is implemented with `Rcpp` and `RcppArmadillo`; the R layer
provides validation, S3 methods, diagnostics, plotting, and documentation.

## Installation

During development, install from the package directory:

```r
install.packages(c("Rcpp", "RcppArmadillo"))
devtools::install("AGCA4extremes")
```

After public release:

```r
install.packages("AGCA4extremes")
```

## Basic Use

```r
library(AGCA4extremes)

data(agca_10d_simulation)

x <- agca_10d_simulation[paste0("X", 1:10)]
fit <- agca(x, k = 250, p = 5)

fit
summary(fit)
plot(fit, type = "variation")
```

The default marginal transformation is rank-Pareto. Larger observations in
each margin are treated as more extreme, so financial return data should be
converted to losses before calling `agca()`.

## Example Diagnostics

The bundled `agca_10d_simulation` data set is generated from a 10-dimensional
heavy-tailed design. Variables `X1`--`X8` share a low-dimensional extremal
mechanism, while `X9` and `X10` contain independent Pareto sources that create
near-axis extreme regimes.

```r
data(agca_10d_simulation)

x <- agca_10d_simulation[paste0("X", 1:10)]
fit <- agca(x, k = 500, p = 4, seed = 1)

agca_rank_summary(fit)
```

Explained variation:

```r
plot(fit, type = "variation")
```

Scores for the first two anchored geodesic components:

```r
cols <- c(shared_low_rank = "#1B9E77", axis_9 = "#5B3A29", axis_10 = "#7570B3")
plot(
  fit$scores[, 1], fit$scores[, 2],
  col = cols[agca_10d_simulation$regime[fit$tail$index]],
  pch = 16,
  xlab = "AGC1 score",
  ylab = "AGC2 score"
)
legend("topright", legend = names(cols), col = cols, pch = 16, bty = "n")
```

Loadings:

```r
plot(fit, type = "loadings", p = 1)
plot(fit, type = "loadings", p = 2)
```

Threshold and anchor diagnostics:

```r
threshold_stability(x, k = c(250, 350, 500, 750), p = 4)
anchor_sensitivity(x, k = 500, p = 4)
```

Bootstrap uncertainty for rank summaries:

```r
boot <- bootstrap_agca(fit, B = 99, ranks = c(1, 2, 4), seed = 1)
summary(boot)
plot(boot, statistic = "variation_explained")
```

## Main Functions

- `agca()` fits the full workflow from data to an AGCA object.
- `rank_pareto()` and `pareto_from_cdf()` perform marginal standardization.
- `tail_directions()` extracts large-radius angular observations.
- `decluster_runs()` performs simple runs declustering for radial extremes.
- `agca_fit_directions()` fits AGCA directly to angular directions.
- `agca_reconstruct()` maps leading anchored coordinates back to the sphere.
- `agca_rank_summary()` reports residual risk and anchored variation explained.
- `bootstrap_agca()` resamples angular directions to quantify sampling
  uncertainty.
- `threshold_stability()` and `anchor_sensitivity()` provide diagnostics.

## Paper Replication

The CRAN package is intentionally lean. Large raw data, generated results,
figures, and full paper workflows should live in a separate
`replicateAGCApaper` repository that depends on `AGCA4extremes`.

## Development Assistance

OpenAI Codex was used as a programming assistant during development, mainly for
code scaffolding, refactoring, documentation, and tests. All methodological
choices, validation, final code, and responsibility for the package remain with
the author.

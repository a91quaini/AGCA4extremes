# AGCA4extremes

<!-- badges: start -->
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

data(ff_portfolio_losses)

x <- ff_portfolio_losses[, -1]
fit <- agca(x, k = 250, p = 5)

fit
summary(fit)
plot(fit, type = "variation")
```

The default marginal transformation is rank-Pareto. Larger observations in
each margin are treated as more extreme, so financial return data should be
converted to losses before calling `agca()`.

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

# Build package-owned simulated example data.
# This script is excluded from the CRAN source tarball by .Rbuildignore.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required to regenerate package data.", call. = FALSE)
}

devtools::load_all(".")

agca_10d_simulation <- simulate_agca_10d(
  n = 10000L,
  seed = 20260627L,
  theta = 0.45,
  tau = 0.25,
  axis9_scale = 1,
  axis10_scale = 1
)

if (!dir.exists("data")) {
  dir.create("data")
}

save(
  agca_10d_simulation,
  file = file.path("data", "agca_10d_simulation.rda"),
  compress = "xz"
)

# Build small package datasets from the paper replication files.
# This script is intentionally excluded from the CRAN source tarball.

source_file <- file.path(
  "..",
  "R",
  "empirics",
  "data_output",
  "sorts_2x3_daily",
  "ff_2x3_sorts_daily.rds"
)

if (!file.exists(source_file)) {
  stop(
    "Cannot find the transformed FF data at ", source_file, ". ",
    "Run this script from the AGCA4extremes package directory inside the ",
    "research repository, or adapt source_file to the replication repo.",
    call. = FALSE
  )
}

ff_source <- readRDS(source_file)
ff_portfolio_losses <- ff_source$complete_losses

if (!dir.exists("data")) {
  dir.create("data")
}

save(
  ff_portfolio_losses,
  file = file.path("data", "ff_portfolio_losses.rda"),
  compress = "xz"
)

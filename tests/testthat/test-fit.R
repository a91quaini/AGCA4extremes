test_that("agca fits and reconstructs angular directions", {
  set.seed(1)
  x <- matrix(stats::rexp(300), ncol = 3)

  fit <- agca(x, k = 30, p = 2)
  expect_s3_class(fit, "agca_fit")
  expect_equal(fit$n_extreme, 30)
  expect_equal(fit$dimension, 3)
  expect_equal(length(fit$eigenvalues), 2)

  ave <- agca_variation_explained(fit)
  expect_equal(length(ave), 2)
  expect_true(all(ave >= 0 & ave <= 1))

  recon <- agca_reconstruct(fit, p = 1)
  expect_equal(dim(recon), dim(fit$g))
  expect_equal(as.numeric(row_norms(recon)), rep(1, nrow(recon)), tolerance = 1e-10)
})

test_that("agca_fit_directions accepts noncanonical anchors", {
  set.seed(2)
  g <- normalize_rows(matrix(stats::rexp(120), ncol = 4))

  fit_f <- agca_fit_directions(g, anchor = "frechet", p = 2, normalize = FALSE)
  fit_p <- agca_fit_directions(g, anchor = "principal", p = 2, normalize = FALSE)

  expect_s3_class(fit_f, "agca_fit")
  expect_s3_class(fit_p, "agca_fit")
  expect_equal(length(fit_f$mu), 4)
  expect_equal(length(fit_p$mu), 4)
})

test_that("declustering returns no more directions than raw exceedances", {
  x <- matrix(seq_len(60), ncol = 3)
  raw <- tail_directions(x, k = 9)
  dec <- decluster_runs(x, k = 9, run = 1)

  expect_lte(nrow(dec$directions), nrow(raw$directions))
  expect_true(length(dec$clusters) >= 1L)
})

test_that("bootstrap and diagnostics return data frames", {
  set.seed(3)
  x <- matrix(stats::rexp(300), ncol = 3)
  fit <- agca(x, k = 30, p = 2)

  boot <- bootstrap_agca(fit, B = 5, seed = 4)
  expect_s3_class(boot, "agca_bootstrap")
  expect_s3_class(summary(boot), "data.frame")

  stab <- threshold_stability(x, k = c(20, 25), p = 2)
  expect_s3_class(stab, "data.frame")
  expect_true(all(c("rank", "k") %in% names(stab)))
})

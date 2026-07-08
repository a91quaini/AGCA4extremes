test_that("simulate_agca_10d returns the documented shape", {
  x <- simulate_agca_10d(n = 50, seed = 1)

  expect_s3_class(x, "data.frame")
  expect_equal(nrow(x), 50)
  expect_true(all(paste0("X", 1:10) %in% names(x)))
  expect_true(is.factor(x$regime))
  expect_true(all(x[paste0("X", 1:10)] > 0))
})

test_that("bundled 10D simulation data are usable by agca", {
  data(agca_10d_simulation)
  x <- agca_10d_simulation[paste0("X", 1:10)]
  fit <- agca(x, k = 100, p = 3)

  expect_s3_class(fit, "agca_fit")
  expect_equal(fit$dimension, 10)
  expect_equal(fit$n_extreme, 100)
})

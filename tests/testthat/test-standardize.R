test_that("rank_pareto returns finite Pareto scores", {
  x <- matrix(c(1, 2, 3, 3, 2, 1), ncol = 2)
  z <- rank_pareto(x)

  expect_equal(dim(z), dim(x))
  expect_true(all(is.finite(z)))
  expect_true(all(z >= 1))
  expect_gt(z[3, 1], z[1, 1])
})

test_that("pareto_from_cdf accepts one CDF or a list of CDFs", {
  x <- matrix(c(0.1, 0.2, 0.3, 0.4), ncol = 2)
  z1 <- pareto_from_cdf(x, cdf = stats::pnorm)
  z2 <- pareto_from_cdf(x, cdf = list(stats::pnorm, stats::pnorm))

  expect_equal(z1, z2)
  expect_true(all(z1 > 1))
})

test_that("agca_standardize supports none", {
  x <- matrix(1:6, ncol = 2)
  expect_equal(agca_standardize(x, margin = "none"), x)
})

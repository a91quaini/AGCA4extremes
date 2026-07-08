test_that("geometry helpers normalize and measure distances", {
  x <- matrix(c(3, 4, 0, 0, 5, 12), ncol = 3, byrow = TRUE)
  g <- normalize_rows(x)

  expect_equal(as.numeric(row_norms(g)), c(1, 1), tolerance = 1e-12)
  expect_equal(sphere_distance(g, g), c(0, 0), tolerance = 1e-12)
})

test_that("anchored departures reject directions outside the open hemisphere", {
  g <- rbind(c(1, 0), c(-1, 0))
  expect_error(anchored_departures(g, c(1, 0)), "open hemisphere")
})

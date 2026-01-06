test_that("skewness and kurtosis return the correct values", {
  expect_equal(skewness(rnorm(100000)), 0, tolerance = 0.05)
  expect_equal(kurtosis(rnorm(100000)), 3, tolerance = 0.05)
})

test_that("Grubbs test and corresponding critical values are correctly computed", {
  y <- replicate(10000, {x <- rnorm(100); grubbs(x)})
  expect_equal(mean(y > crit.grubbs(100, 0.05)), 0.05, tolerance = 0.05)
})
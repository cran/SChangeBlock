test_that("GMD correctly computes Gini's mean difference", {
  x <- rnorm(10)
  
  expect_equal(GMD(x), mean(dist(x)))
})

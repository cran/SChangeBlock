test_that("autocov has the correct format", {
  X <- genField(c(20, 30))
  Sigma <- autocov(X, b = c(4, 4))
  Sigma1 <- autocov(X, b = 4, direction = 1)
  Sigma2 <- autocov(X, b = 4, direction = 2)
  
  expect_equal(dim(Sigma), c(600, 600))
  expect_equal(dim(Sigma1), c(30, 30))
  expect_equal(dim(Sigma2), c(20, 20))
  
  expect_equal(Sigma[1, 6], 0)
  expect_equal(Sigma[20, 1], 0)
  expect_true(Sigma[21, 1] != 0)
})
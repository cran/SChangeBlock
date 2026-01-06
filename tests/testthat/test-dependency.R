test_that("dependency matrix Phi is correct", {
  thetaAR <- genTheta(1, 0.4, structure = "AR")
  
  expect_true("matrix" %in% class(thetaAR))
  expect_equal(dim(thetaAR), c(3, 3))
  theta_test <- matrix(c(0.4^sqrt(2), 0.4, 0.4^sqrt(2), 
                         0.4, 1, 0.4, 
                         0.4^sqrt(2), 0.4, 0.4^sqrt(2)), ncol = 3)
  theta_test <- theta_test / sqrt(sum(theta_test^2))
  expect_equal(thetaAR, theta_test)
  
  thetaMA <- genTheta(1, 0.4, structure = "MA")
  theta_test <- matrix(c(0.4^2, 0.4, 0.4^2, 
                         0.4, 1, 0.4, 
                         0.4^2, 0.4, 0.4^2), ncol = 3)
  expect_equal(thetaMA, theta_test)
  
  phi <- genTheta(2, c(0.8, 0.6, 0.4, 0.2))
  expect_equal(dim(phi), c(5, 5))
  expect_equal(phi, matrix(c(0.2, 0.4, 0.6, 0.4, 0.2, 
                             0.4, 0.6, 0.8, 0.6, 0.4, 
                             0.6, 0.8, 1, 0.8, 0.6, 
                             0.4, 0.6, 0.8, 0.6, 0.4, 
                             0.2, 0.4, 0.6, 0.4, 0.2), ncol = 5))
})

test_that("dependency is correctly added to the data", 
{
  E <- matrix(rnorm(25), ncol = 5)
  theta <- genTheta(2, 0.4)
  X1 <- dependency(E, q_ = 2, param_ = 0.4)
  X2 <- dependency(E, Theta_ = theta)
  Xtest <- sum(E * theta)
  
  expect_equal(X1, X2)
  expect_true("matrix" %in% class(X1))
  expect_equal(dim(X1), c(1, 1))
  expect_equal(X1[1], Xtest)
  expect_equal(X2[1], Xtest)
})
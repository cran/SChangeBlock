test_that("bandwidth estimation returns the correct format", {
  x <- arima.sim(list(ar = 0.5), 100)
  X <- genField(c(50, 50), Theta = genTheta(10, 0.4))

  expect_equal(length(bandwidth(x, 0.3, 0.3)), 1)
  expect_equal(length(bandwidth(X, 0.3, 0.3)), 2)
})


test_that("de-correlation is working correctly", {
  X <- genField(c(20, 20), Theta = genTheta(10, 0.3))
  Y <- decorr(X, lags = bandwidth(X, 0.3, 0.3))
  
  expect_equal(class(X), class(Y))
  
  x <- arima.sim(model = list(ar = 0.5), 100)
  y <- decorr(x, bandwidth(x, 0.3, 0.3))
  
  expect_equal(class(x), class(y))
  
  ## check if arima parameter is significant (to 10% level)
  arima.x <- arima(x, c(1, 0, 0))
  arima.y <- arima(y, c(1, 0, 0))
  
  expect_lt((1 - pnorm(abs(arima.x$coef[1]) / sqrt(arima.x$var.coef[1])))*2, 0.1)
  expect_gt((1 - pnorm(abs(arima.y$coef[1]) / sqrt(arima.y$var.coef[1])))*2, 0.1)
})


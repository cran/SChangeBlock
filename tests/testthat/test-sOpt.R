test_that("sOpt returns the correct exponent", {
  n <- c(10, 20, 50, 100, 200, 400)
  s.opt <- sOpt(n, 0.6)
  l <- round(n^s.opt)
  b <- floor(n / l)
  
  expect_equal(l * b, n)
})
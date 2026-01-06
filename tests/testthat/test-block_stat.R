test_that("block_stat produces the correct format", {
  X <- genField(c(50, 50))
  
  y1 <- block_stat(X, 0.6, "gmd")
  y2 <- block_stat(X, 0.6, "var")
  y3 <- block_stat(X, 0.6, "jb")
  y4 <- block_stat(X, 0.6, "ks")
  y5 <- block_stat(X, 0.6, "grubbs")
  y6 <- block_stat(X, 0.6, "ANOVA")
  
  expect_equal(length(y1), 1)
  expect_equal(length(y2), 1)
  expect_equal(length(y3), 1)
  expect_equal(length(y4), 1)
  expect_equal(length(y5), 1)
  expect_equal(length(y6), 1)
  
  expect_equal(attributes(y5), list(n = 25, blocksize = c(10, 10)))
  expect_equal(attributes(y6), list(k = 25, N = 50^2, blocksize = c(10, 10)))
})

test_that("block_stat and block_pValue fit together", {
  tn <- replicate(1000, 
  {
    X <- genField(c(50, 50))
    block_stat(X, 0.6, "gmd")
  })
  expect_lt(mean(block_pValue(tn, "gmd") < 0.05), 0.1)
  
  tn <- replicate(1000, 
  {
    X <- genField(c(50, 50))
    block_stat(X, 0.6, "var")
  })
  expect_lt(mean(block_pValue(tn, "var") < 0.05), 0.1)
  
})
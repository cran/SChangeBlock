test_that("block_test produces the correct output",
{
  x <- rnorm(c(50, 50))
  y <- block_test(x, s = 0.6, fun = "gmd")
  
  expect_true(is.list(y))
  expect_equal(class(y), "htest")
  
  expect_true(all(c("alternative", "method", "data.name",
                    "statistic", "p.value") %in% names(y)))
})

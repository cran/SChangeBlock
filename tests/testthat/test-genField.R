test_that("changeRegion returns the correct format", {
  expect_equal(class(changeRegion(c(50, 50), type = "1a")), "integer")
  expect_equal(class(changeRegion(c(50, 50), type = "1b")), "integer")
  expect_equal(class(changeRegion(c(50, 50), type = "1c")), "integer")
  expect_equal(class(changeRegion(c(50, 50), type = 2)), "integer")
  expect_equal(class(changeRegion(c(50, 50), type = 4)), "integer")
  expect_equal(class(changeRegion(c(50, 50), type = 5)), "integer")
  
  expect_equal(class(changeRegion(c(50, 50), type = 3)), c("matrix", "array"))
  expect_length(changeRegion(c(50, 50), type = 3), 50^2)
})

test_that("genField returns the correct format", {
  X <- genField(c(50, 40))
  
  expect_equal(dim(X), c(50, 40))
  expect_equal(class(X), "RandomField")
  
  x <- genField(50)
  
  expect_equal(length(x), 50)
  expect_equal(class(x), "RandomField")
})

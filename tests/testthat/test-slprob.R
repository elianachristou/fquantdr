test_that("the function returns an error when y is not a vector", {
  y <- -10000
  # This statement should throw an error
  H <- 3
  expect_no_error(slprob(y, H))
})

test_that("H has to be a positive integer", {
  y <- c(4, 7, 8, 9, 17)
  H <- -3
  expect_error(slprob(y, H))
})

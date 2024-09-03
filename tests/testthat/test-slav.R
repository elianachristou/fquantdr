library(fquantdr)

test_that("The input is a matrix", {
  x <- 3
  y <- c(2, 4, 5, 7, 8, 9)
  H <- 4
  expect_error(slav(x, y, H))
})

test_that("the input for y is a vector", {
  y <- c(2, 5, 1, 6, 9, 6, 0)
  H <- 3
  expect_no_error(discretize(y, H))
})

test_that("the input for H positive whole number", {
  # This contains the same issue as slprob, the function allows scalar values
  y <- c(2, 4, 5, 2, 8, 1)
  H <- 0
  expect_error(discretize(y, H))
})
# ERROR this function returns the matrix in an order of "how the y vector is classified
# by discretize. Instead of returning the matrix it instead returns the matrix in
# a different order

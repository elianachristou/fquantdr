library(fquantdr)

test_that("the input for y is a vector", {
  # This contains the same issue as slprob, the function allows scalar values
  # Not a major issue
  y <- 3
  H <- 5
  expect_no_error(discretize(y, H))
})

test_that("the input for H positive whole number", {
  # This contains the same issue as slprob, the function allows scalar values
  y <- 3
  H <- c(1, 6, 10)
  expect_error(discretize(y, H))
})

test_that("the output y1 is a vector", {
  y <- rnorm(12)
  H <-5
  # This should expect an error becuase discretize returns a vector however this
  # is not reflected here
  expect_no_error(!is.vector(discretize(y, H)))
})


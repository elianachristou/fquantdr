library(fquantdr)

test_that("H is a number less than the length of y", {
  y <- rnorm(4)
  H <- 5
  expect_error(discretize(y, H))
})

test_that("H is a single number", {
  y <- c(1, 3, 5, 9, 4)
  H <- c(1, 6, 10)
  expect_error(discretize(y, H))
})

test_that("H is a positive integer number", {
  y <- c(1, 3, 5, 9, 4)
  H <- -2.2
  expect_error(discretize(y, H))
})



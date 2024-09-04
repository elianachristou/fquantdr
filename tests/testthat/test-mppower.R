library(fquantdr)

test_that("input a is a square matrix", {
  a <- matrix(rnorm(12), nrow = 2, ncol = 6)
  alpha <- 2
  expect_error(mppower(a, alpha))
})

test_that("alpha is numeric", {
  a <- matrix(rnorm(36), nrow = 6, ncol = 6)
  alpha <- "how"
  expect_error(mppower(a, alpha))
})

test_that("alpha is a single number", {
  a <- matrix(rnorm(36), nrow = 6, ncol = 6)
  alpha <- c(-1, 2)
  expect_error(mppower(a, alpha))
})

test_that("epsilon is a real number", {
  a <- matrix(rnorm(36), nrow = 6, ncol = 6)
  alpha <- 3
  epsilon <- 3 + 2i
  expect_error(mppower(a, alpha, epsilon))
})

test_that("epsilon is a nonnegative real number", {
  a <- matrix(rnorm(36), nrow = 6, ncol = 6)
  alpha <- 3
  epsilon <- -3
  expect_error(mppower(a, alpha, epsilon))
})

#test_that("mppower gives the same output as matpower if epsion is 0", {
#  a <- matrix(c(6, 4, 8, 2, 5, 9, 3, 1, 7), nrow = 3, ncol = 3)
#  alpha <- 2
#  epsilon <- 0
#  expect_equal(mppower(a, alpha, epsilon), matpower(a, alpha))
#})

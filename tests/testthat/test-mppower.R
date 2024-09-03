library(fquantdr)

test_that("input a is a square matrix", {
  a <- matrix(rnorm(12), nrow = 2, ncol = 6)
  alpha <- 2
  expect_error(mppower(a, alpha))
})

test_that("alpha is a positive scalar value", {
  a <- matrix(rnorm(36), nrow = 6, ncol = 6)
  alpha <- "how"
  expect_error(mppower(a, alpha))
})

test_that("epsilon inputs are not complex", {
  a <- matrix(rnorm(36), nrow = 6, ncol = 6)
  alpha <- 3
  epsilon <- 3 + 2i
  expect_error(mppower(a, alpha, epsilon))
})

test_that("epsilon inputs are not positive", {
  a <- matrix(rnorm(36), nrow = 6, ncol = 6)
  alpha <- 3
  epsilon <- -3
  expect_error(mppower(a, alpha, epsilon))
})

test_that("ignore must be a numeric value", {
  a <- matrix(rnorm(36), nrow = 6, ncol = 6)
  alpha <- 3
  epsilon <- 3
  ignore <- 3 + 2i
  expect_error(mppower(a, alpha, epsilon, ignore))
})

test_that("output a is a square matrix", {
  a <- matrix(rnorm(12), nrow = 2, ncol = 6)
  alpha <- 2
  # How do I test if the output is a matrix?
  # This does not return an error however since I said expect error and it
  # Does not create an error it should return an error
  expect_error(is.matrix(mppower(a, alpha)))
})

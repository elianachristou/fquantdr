library(fquantdr)

test_that("y is a univariate response", {
  n <- 100
  nbasis <- 10
  p <- 3
  basis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
  coef_matrix <- array(rnorm(n * nbasis * p), dim = c(nbasis, n, p))
  xfd <- fda::fd(coef_matrix, basis)
  y <- matrix(1:6, nrow = 6, ncol = 1)
  expect_error(sonf(y, xfd))
})

test_that("y is of length n", {
  n <- 100
  m <- 80
  nbasis <- 8
  p <- 5
  basis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
  coef_matrix <- array(rnorm(n * nbasis * p), dim = c(nbasis, n, p))
  xfd <- fda::fd(coef_matrix, basis)
  y <- rnorm(m)
  expect_error(sonf(y, xfd))
})

test_that("xfd is a functional data object of the class 'fd'", {
  n <- 100
  nbasis <- 12
  p <- 3
  basis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
  coef_matrix <- array(rnorm(n * nbasis * p), dim = c(nbasis, n, p))
  y <- rnorm(n)
  expect_error(sonf(y, coef_matrix))
})

test_that("xfd is a 3-dimensional array", {
  n <- 100
  nbasis <- 7
  p <- 2
  basis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
  coef_matrix <- array(rnorm(n * nbasis), dim = c(nbasis, n))
  xfd <- fda::fd(coef_matrix, basis)
  y <- rnorm(n)
  expect_error(sonf(y, xfd))
})


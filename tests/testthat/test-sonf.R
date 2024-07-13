test_that("y is a univariate response", {
  set.seed(1)
  n <- 100
  nbasis <- 10
  p <- 3
  basis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
  coef_matrix <- array(rnorm(n * nbasis * p), dim = c(nbasis, n, p))
  fdobj <- fda::fd(coef_matrix, basis)
  # Creates a multivariate vector
  y <- matrix(1:6, nrow = 6, ncol = 1)
  expect_error(sonf(y, fdobj))
})

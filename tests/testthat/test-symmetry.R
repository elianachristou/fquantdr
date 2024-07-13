library(fquantdr)

test_that("The input is a matrix", {
  sym.mat <- 3
  expect_error(symmetry(sym.mat))
})

test_that("The matrix is square", {
  sym.mat <- matrix(rnorm(6), nrow = 2, ncol = 3)
  expect_error(symmetry(sym.mat))
})

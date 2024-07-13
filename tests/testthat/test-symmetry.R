library(fquantdr)

test_that("The input is a matrix", {
  sym.mat <- 3
  expect_error(symmetry(sym.mat))
})

library(fquantdr)

test_that("n is a positive integer", {
  n <- -2.1
  expect_error(qmat(n))
})


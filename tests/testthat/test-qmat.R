library(fquantdr)

test_that("n is a positive integer", {
  n <- -2.1
  expect_error(qmat(n))
})

test_that("n can be a floating value",{
  n <- 1.00
  expect_no_error(qmat(n))
})

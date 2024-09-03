library(fquantdr)

test_that("x is a matrix", {
  x <- 3
  y <- c(2, 4, 5, 7, 8, 9)
  H <- 4
  ydis <- discretize(y, H)
  expect_error(slav(x, ydis, H))
})


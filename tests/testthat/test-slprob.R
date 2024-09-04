library(fquantdr)

test_that("H is a positive integer", {
  ydis <- c(4, 7, 8, 9, 17)
  H <- -3
  expect_error(slprob(ydis, H))
})

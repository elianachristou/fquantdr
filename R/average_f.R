average_f <- function(x) {
  # this finds the average of x
  output1 <- sum(x) / length(x)
  # this finds the sum of the squared values divide by n
  output2 <- sum(x^2) / length(x)
  # this calculates the variance
  output3 <- sum(x - output1)^2 / (length(x) - 1)
  list(output1 = output1, output2 = output2, output3 = output3)
}

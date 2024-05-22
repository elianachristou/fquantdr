average_f <- function(x) {
  # this finds the average of x
  output1 <- sum(x) / length(x)
  # this finds the sum of the squared values divide by n
  output2 <- sum(x^2) / length(x)
  list(output1 = output1, output2 = output2)
}

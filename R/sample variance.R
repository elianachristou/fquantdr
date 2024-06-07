# Function to calculate sample variance
calculate_sample_variance <- function(x) {
  n <- length(x)
  if (n <= 1) {
    stop("Sample size must be greater than 1")
  }
  mean_x <- mean(x)
  variance <- sum((x - mean_x)^2) / (n - 1)
  return(variance)
}

calculate_median <- function(x) {
  if (!is.numeric(x)) {
    stop("Input must be a numeric vector.")
  }
  return(median(x))
}

# Example usage:
example_vector <- c(1, 2, 3, 4, 5)
A <- calculate_median(example_vector)
A

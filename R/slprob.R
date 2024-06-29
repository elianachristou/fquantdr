#' Slice Probabilities
#'
#' \code{slprob} calculates the probability that a given response will be in
#' each slice.
#'
#' This function takes in a discretized y vector and a vector determining the
#' values at each slice (usually a range 1:H where H is the number of slices)
#' to determine the proportion of the observed responses in each slice to
#' determine the probablity that a given y will be within that slice.
#'
#' @param y The discretized vector for the response variable
#' @param H The number of slices used
#'
#' @return A vector of probabilities that a response will be in each slice
#'
#' @noRd
#' @examples
#' slprob(c(1,1,2,3,3,2,1,2), 3)
#'
slprob <- function(y, H) {
  if (!is.vector(y)) {
    stop("y must be a vector.")
  }
  if (!is.numeric(H) | H <= 0 | as.integer(H) != H) {
    stop("H must be a positive integer.")
  }
  for (i in y) {
    if(!(i %in% 1:H)) {
      stop(paste("y must be discretized such that all values in y are in",
                 "yunit. The value", i, "is not in yunit."))
    }
  }

  # Get length of y vector
  n <- length(y)
  yunit <- 1:H
  # Define new vector for output
  out <- rep(0, H)
  for (i in yunit) {
    # Probabilty = number of points in slice / number of points
    out[i] <- length(y[y == yunit[i]]) / n
  }
  # Return output vector
  out
}

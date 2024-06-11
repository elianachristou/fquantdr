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
#' @param yunit A vector defining the slices used
#'
#' @return A vector of probabilities that a response will be in each slice
#'
#' @noRd
#' @examples
#' slprob(c(1,1,2,3,3,2,1,2),1:3)
slprob <- function(y, yunit) {
  # Get length of y vector
  n <- length(y)
  # Get number of slices from yunit
  nslice <- length(yunit)
  # Define new vector for output
  out <- rep(0, nslice)
  for (i in 1:nslice) {
    # Probabilty = number of points in slice / number of points
    out[i] <- length(y[y == yunit[i]]) / n
  }
  # Return output vector
  out
}

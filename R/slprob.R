#' Slice Probabilities
#'
#' \code{slprob} calculates the proportions of observations falling within
#' each slice based on a discretized response vector.
#'
#' This function calculates the slice probabilities for a given discrete
#' numeric vector \code{ydis}, where each value represents the slice
#' assignment of an observation.  The total number of slices is specified by
#' \code{H}.  The function returns an \code{H}-dimensional vector containing
#' the proportion of observations in each slice.
#'
#' @param ydis A discrete numeric vector that contains labels that define the
#' slice each value of a vector `y` is in. This can be obtained using the
#' `discretize` function.
#' @param H The number of slices.
#'
#' @return A `H`-dimensional vector of proportions of observations that fall
#' within each slice.
#'
#' @seealso \code{\link{discretize}}, \code{\link{slav}}
#' @noRd
#' @examples
#' y <- rnorm(100)
#' H <- 3
#' ydis <- discretize(y, H)
#' slprob(ydis, H)
#'
slprob <- function(ydis, H) {

  # Check if ydis is a vector
  if (!is.vector(ydis)) {
    stop("ydis must be a vector.")
  }

  # Check if H is a positive integer
  if (H <= 0 | H != floor(H)) {
    stop("H must be a positive integer.")
  }

  # Check if H is one number
  if (length(H) > 1) {
    stop("H must be one number.")
  }

  # Check if H is greater than 1 and less than the length of ydis
  if (H == 1 | H >= length(ydis)) {
    stop("H must be an integer that is at least 2 and less than the
         length of ydis.")
  }

  # Define the parameters
  n <- length(ydis)
  yunit <- 1:H

  # Calculate the proportions
  out <- rep(0, H)
  for (i in yunit) {
    # Probabilty = number of points in slice / number of points
    out[i] <- length(ydis[ydis == yunit[i]]) / n
  }
  out
}

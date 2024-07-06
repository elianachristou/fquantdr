#' Slice Probabilities
#'
#' \code{slprob} calculates the proportion of observations that fall within
#' each slice, where the number of slices is given by H.
#'
#' @param y A numeric vector.
#' @param H The number of slices.
#'
#' @return A vector of proportion of observations that fall within each slice.
#'
#'@seealso \code{\link{discretize}}, \code{\link{slav}}
#' @noRd
#' @examples
#' y <- rnorm(100)
#' H <- 3
#' slprob(y, H)
#'
slprob <- function(y, H) {

  # Check if y is a vector
  if (!is.vector(y)) {
    stop("y must be a vector.")
  }

  # Check if H is a positive integer
  if (H <= 0 | H != floor(H)) {
    stop("H must be a positive integer.")
  }

  # Define the parameters
  n <- length(y)
  yunit <- 1:H
  ydis <- discretize(y, H)

  # Calculate the proportions
  out <- rep(0, H)
  for (i in yunit) {
    # Probabilty = number of points in slice / number of points
    out[i] <- length(ydis[ydis == yunit[i]]) / n
  }
  out
}

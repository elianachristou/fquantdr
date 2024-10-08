#' Discretize conversion of a continuous numeric vector
#'
#' \code{discretize} converts a continuous numeric vector into a discrete
#' vector with labels corresponding to different slices of the data.
#'
#' This function converts a continuous numeric vector into a discrete vector
#' by slicing the data and assigning a label to each value according to which
#' slice it is in.  The idea stems from Li (1991), who proposed sliced inverse
#' regression (SIR), a dimension reduction technique.  This function is
#' necessary for this package when performing the extension of SIR to
#' functional predictors, introduced by Ferr&#233; and Yao (2003).
#'
#' @param y A (continuous or discrete) vector.  If y is discrete, a small
#'     amount of noise is added to make it continuous.
#' @param H The number of slices.
#'
#' @return A vector defining the slices each value of y corresponds to.
#'
#' @references Li, K.-C. (1991) Sliced Inverse Regression for Dimension
#' Reduction.  \emph{Journal of the American Statistical Association},
#' 86(414), 316-327.
#'
#' Ferr&#233;, L, and Yao, F. (2003) Function Sliced Inverse Regression Analysis.
#' \emph{Statistics}, 37(6), 475-488.
#'
#' @noRd
#' @examples
#' # Example 1
#' y <-  c(2.5, 3.6, 1.2, 4.8, 2.9)
#' H <- 3
#' discretize(y, H)
#'
#' # Example 2
#' y <- rnorm(100)
#' H <- 4
#' discretize(y, H)
#'
discretize <- function(y, H) {

  # Check if y is a vector
  if (!is.vector(y)) {
    stop("y must be a vector.")
  }

  # Check if H is one number
  if (length(H) > 1) {
    stop("H must be one number.")
  }

  # Check if H is a positive integer
  if (H != round(H) | H <= 0) {
    stop("H must be a positive integer.")
  }

  # Check if H is greater than 1 and less than the length of y
  if (H == 1 | H >= length(y)) {
    stop("H must be an integer that is at least 2 and less than the
         length of y.")
  }

  # define the parameters
  n <- length(y)
  yunit <- 1:H

  # Add small amount of noise to y
  y <- y + .00001 * mean(y) * stats::rnorm(n)

  # Order y values in ascending order
  yord <- y[order(y)]
  # Calculate the approximate number of data points per slice
  nwidth <- floor(n / H)

  # Set each division point to the values of yord that represent
  # the boundaries between slices
  divpt <- rep(0, H - 1)
  for(i in 1:(H - 1)) {
    divpt[i] <- yord[i * nwidth + 1]
  }

  y1 <- rep(0, n)
  # Assign slice labels to the upper boundary slice
  y1[y >= divpt[H - 1]] <- H
  # Assign slice labels to the lower boundary slice
  y1[y < divpt[1]] <- 1
  # Assign slices labels to intermediate slices
  for(i in 2:(H - 1)) {
    y1[(y >= divpt[i - 1]) & (y < divpt[i])] <- i
  }

  # Return discretized y
  y1
}

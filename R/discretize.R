#' Discretize conversion of a continuous numeric vector
#'
#' \code{discretize} partitions a continuous numeric vector into a discrete
#' vector with labels corresponding to different slices of the data.
#'
#' This function discretizes a continuous numeric vector into a discrete vector
#' by dividing its rage into \code{H} approximately equal-sized slices and
#' labeling each observation based on its slice assignment.  If \code{y} is
#' discrete, a small amount of noise is added to prevent ties. The idea
#' originates from Sliced Inverse Regression (SIR) proposed by Li (1991) and is
#' used in functional SIR, as extended by Ferr&#233; and Yao (2003).
#'
#' @param y A (continuous or discrete) numeric vector.  If y is discrete, a small
#'     amount of noise is added to ensure continuity.
#' @param H A positive integer specifying the number of slices.
#'
#' @return A numeric vector of the same length as \code{y}, with integer labels
#'     \code{1} to \code{H}, indicating the slice assignment of each value.
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
#' discretized_y <- discretize(y, H)
#' table(discretized_y)
#'
discretize <- function(y, H) {

  # Check if y is a vector
  if (!is.vector(y) || !is.numeric(y)) {
    stop("y must be a numeric vector.")
  }

  # Check if H is a single number
  if (length(H) > 1) {
    stop("H must be a single number.")
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

  # Add small noise if y is discrete (prevents ties)
  if (length(unique(y)) < length(y) * 0.75) {
  y <- y + .00001 * mean(y) * stats::rnorm(n)
  }

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

#' Discretize conversion of a continuous numeric vector
#'
#' \code{discretize} converts a continuous numeric vector into a discrete
#' vector with labels corresponding to different slices of the data.
#'
#' This function converts a continuous numeric vector into a discrete vector
#' by slicing the data and assigning a label to each value according to which
#' slice it is in.  The idea stems from Li (1991), who proposed sliced inverse
#' regression (SIR), a dimension reduction technique.  This function is
#' necessary when performing the extension of SIR to functional predictors,
#' introduced by Ferr\'e and Yao (2003).
#'
#' @param y A (continuous or discrete) vector
#' @param H The number of slices
#'
#' @return A vector defining the slices each value of y corresponds to
#'
#' @references Li, K.-C. (1991) Sliced Inverse Regression for Dimension Reduction.
#'     \emph{Journal of the American Statistical Association}, 86(414), 316-327.
#'    Ferr\'e, L, and Yao, F. (2003) Function Sliced Inverse Regression
#'    Analysis. \emph{Statistics}, 37(6), 475-488.
#' @noRd
#' @examples
#' y <-  c(1, 2, 3, 2)
#' H <- 3
#' discretize(y, H)
#'
discretize <- function(y, H) {

  # Check if y is a vector
  if (!is.vector(y)) {
    stop("y must be a vector.")
  }

  # Check if H is a number
  if (!is.numeric(H)) {
    stop("H must be a number.")
  }

  n <- length(y)
  yunit <- 1:H

  # Add small amount of noise to y
  y <- y + .00001 * mean(y) * stats::rnorm(n)
  nsli <- length(yunit)

  # Order y values in ascending order
  yord <- y[order(y)]
  n <- length(y)
  nwidth <- floor(n / nsli)

  # Instantiate vector of division points between slices
  divpt <- rep(0, nsli - 1)

  # Set each division point to the values of yord that represent the boundaries
  # between slices
  for(i in 1:(nsli - 1)) {
    divpt[i] <- yord[i * nwidth + 1]
  }

  y1 <- rep(0, n)
  # Assign slice labels to the upper boundary slice
  y1[y >= divpt[nsli - 1]] <- nsli
  # Assign slice labels to the lower boundary slice
  y1[y < divpt[1]] <- 1
  # Assign slices labels to intermediate slices
  for(i in 2:(nsli - 1)) {
    y1[(y >= divpt[i - 1]) & (y < divpt[i])] <- i
  }
  # Return discretized y
  y1
}

#' Compute slice average
#'
#' \code{slav} computes the slice average.
#'
#' Creates a matrix of slice average by using x coordinates, a discretized
#' response, and number of slices. For each slice, the rows of x coefficients
#' where y equals the current slice value of y unit at i get the mean function
#' applied to them column-wise to compute the average of x coefficient for each
#' slice.
#'
#' @param xcoefs Coordinates of X
#' @param y The discretized vector for the response variable
#' @param yunit A vector defining the slices used
#'
#' @return Matrix of slice average
#' @noRd
#'
#' @examples
#' xcoefs <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
#'   nrow = 4, byrow= TRUE)
#' y <- c(1, 2, 3, 2)
#' yunit <- unique(y)
#' slav(xcoefs, y, yunit)
slav <- function(xcoefs, y, yunit) {
  # Get dimensions of xcoefs
  n <- nrow(xcoefs)
  K <- ncol(xcoefs)
  # Get number of slices
  nslice <- length(yunit)
  # Initialize matrix of mean coordinates for each slice
  xcoefsgy <- matrix(0, nslice, K)
  for(i in 1:nslice) {
    # Store predictors for each slice
    predictors <- xcoefs[y == yunit[i], ]
    # If vector, convert to horizontal matrix
    if (is.vector(predictors)) {
      predictors <- t(as.matrix(predictors))
    }
    # Calculate mean for each slice
    xcoefsgy[i, ] <- apply(predictors, 2, mean)
  }
  # Return mean coordinates for the predictors
  xcoefsgy
}

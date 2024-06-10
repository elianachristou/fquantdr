#' Compute slice average
#'
#' \code{slav} computes slice average
#'
#' Creates a matrix of slice average by using x coordinates, a discretized
#' response, and number of slices. For each slice, the rows of x coefficients
#' where y equals the current slice value of y unit at i get the mean function
#' applied to them column-wise to compute the average of x coefficient for each
#' slice.
#'
#' @param xcoefs coordinates of X
#' @param y discretized response
#' @param yunit number of slices
#'
#' @return matrix of slice average
#' @noRd
#'
#' @examples
#' xcoefs <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
#'   nrow = 4, byrow= TRUE)
#' y <- (1, 2 , 3, 2)
#' yunit <- unique(y)
slav <- function(xcoefs, y, yunit) {
  n <- nrow(xcoefs)
  K <- ncol(xcoefs)
  nslice <- length(yunit)
  xcoefsgy <- matrix(0, nslice, K)
  for(i in 1:nslice) {
    xcoefsgy[i, ] <- apply(xcoefs[y==yunit[i], ], 2, mean) }
  return(xcoefsgy)
}

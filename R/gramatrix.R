#' Computes an inner product matrix

#' \code{gramatrix} calculates the inner product of a spline and forms a matrix


#' This function computes the center gram matrix from the inner products of
#' databasis functions. That inner product is taken and used as an input
#' for the qmat function to return a matrix where the sum of each row and column
#' is zero.

#' @param databasis a spline object that is a representation of data as a
#' smoothed function
#' @param K the number of elements in the basis

#' @return A K x K matrix that has been scaled by its inner product

#' @noRd

#' databasis is a spline object



#' @examples
#' K <- 3
#' x <- seq(0, 10, length.out = 5)
#' b_spline <- bs(x, degree = 3, knots = c(3, 7))
#' databasis <- b_spline
#' grammatrix(K, databasis)


gramatrix <- function(K, databasis) {
  # compatability checks for both K and databasis
  if (!is.numeric(K) || K <= 0 || K != as.integer(K)) {
    stop("K must be a positive integer.")
  }
  if (!inherits(databasis, "basisfd")) {
    stop("databasis must be a 'basisfd' object.")
  }

  # Kernmat is assigned a K x K matrix of the inner products
  kernmat <- inprod(databasis, databasis)

  # Checks to make sure inprod returned a matrix
  if (!is.matrix(kernmat)) {
    stop("The inprod function has not returned a matrix.")
  }
  # Computes a K x K matrix using the qmat function for centering
  Gmat <- qmat(K) %*% kernmat %*% qmat(K)

  # returns a centered inner product matrix
  return(Gmat)
}

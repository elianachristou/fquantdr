#' Compute the Gram matrix
#'
#' \code{gramatrix} calculates the inner products of basis functions.
#'
#' This function computes the centered Gram matrix from the inner products
#' of given basis functions, commonly used in functional data analysis.
#'
#' @param K A positive integer representing the number of elements in the
#'     basis, usually termed as `nbasis`.  It must be a positive integer.
#' @param databasis A basis object that is a representation of data as a
#'     smoothed function.  It must be an object of class \code{"basisfd"},
#'     typically created using functions from the `fda` package.
#'
#' @return A centered \code{K x K} inner product matrix, representing the
#'     Gram matrix.
#'
#' @noRd
#' @examples
#' # Example 1: Compute the Gram matrix for a B-spline basis
#' K <- 4
#' databasis <- fda::create.bspline.basis(rangeval = c(0, 10), nbasis = K)
#' gramatrix(K, databasis)
#'
#' \dontrun{
#' # Example 2: Error message when databasis is not a basisfd object
#' K <- 5
#' x <- seq(0, 10, length.out = 5)
#' databasis <- splines::bs(x, degree = 3, knots = c(3, 7))
#' gramatrix(K, databasis)
#' }
#'
gramatrix <- function(K, databasis) {

  # Check if K is a positive integer
  if (K <= 0 | K != round(K)) {
    stop("Parameter 'K' must be a positive integer.")
  }

  # Check if databasis is a basisfd object
  if (!inherits(databasis, "basisfd")) {
    stop("Input 'databasis' must be a 'basisfd' object.")
  }

  # Calculate the K x K matrix of the inner products
  kernmat <- fda::inprod(databasis, databasis)

  # Calculate the centered Gram matrix
  Gmat <- qmat(K) %*% kernmat %*% qmat(K)
  Gmat
}

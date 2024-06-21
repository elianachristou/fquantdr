#' Compute the Gram matrix
#'
#' \code{gramatrix} calculates the inner products of basis functions.
#'
#' This function computes the centered Gram matrix from the inner products
#' of given basis functions, commonly used in functional data analysis.
#'
#' @param K The number of elements in the basis, which must be a positive
#'     integer.
#' @param databasis A basis object that is a representation of data as a
#'     smoothed function.  It must be an object of class \code{"basisfd"},
#'     typically created using functions from the `fda` package.
#'
#' @return A centered \code{K x K} inner product matrix, representing the Gram
#'     matrix.
#'
#' @noRd
#' @examples
#' # Example 1
#' K <- 4
#' x <- seq(0, 10, length.out = 5)
#' databasis <- fda::create.bspline.basis(rangeval = c(0, 10), nbasis = K,
#'     norder = 4)
#' gramatrix(K, databasis)
#'
#' \dontrun{
#' Example 2
#' x <- seq(0, 10, length.out = 5)
#' databasis <- bs(x, degree = 3, knots = c(3, 7))
#' gramatrix(K, databasis)
#' }
#'
gramatrix <- function(K, databasis) {

  # compatibility checks for both K and databasis
  if (!is.numeric(K) | K <= 0 | K != round(K)) {
    stop("K must be a positive integer.")
  }

  if (!inherits(databasis, "basisfd")) {
    stop("databasis must be a 'basisfd' object.")
  }

  # Kernmat is assigned a K x K matrix of the inner products
  kernmat <- fda::inprod(databasis, databasis)

  # Computes a K x K matrix using the qmat function for centering
  Gmat <- qmat(K) %*% kernmat %*% qmat(K)

  # returns a centered inner product matrix
  Gmat
}

#' Gram matrix
#'
#' \code{gramatrix} calculates the inner products of basis functions.
#'
#' This function computes a centered Gram matrix from the inner products
#' of a given basis object, which is often utilized in functional data
#' analysis.
#'
#' @param K The number of elements in the basis, which must be a positive
#'     integer.
#' @param databasis A basis object that is a representation of data as a
#' smoothed function.
#'
#' @return A centered K x K inner product matrix, representing the Gram
#'     matrix.
#'
#' @noRd
#' @examples
#' K <- 3
#' x <- seq(0, 10, length.out = 5)
#' b_spline <- bs(x, degree = 3, knots = c(3, 7))
#' databasis <- b_spline
#' grammatrix(K, databasis)
#'
gramatrix <- function(K, databasis) {
  # compatability checks for both K and databasis
  if (!is.numeric(K) || K <= 0 || K != as.integer(K)) {
    stop("K must be a positive integer.")
  }
  if (!inherits(databasis, "basisfd")) {
    stop("databasis must be a 'basisfd' object.")
  }

  # Kernmat is assigned a K x K matrix of the inner products
  kernmat <- fda::inprod(databasis, databasis)

  # Checks to make sure inprod returned a matrix
  if (!is.matrix(kernmat)) {
    stop("The inprod function has not returned a matrix.")
  }
  # Computes a K x K matrix using the qmat function for centering
  Gmat <- qmat(K) %*% kernmat %*% qmat(K)

  # returns a centered inner product matrix
  return(Gmat)
}

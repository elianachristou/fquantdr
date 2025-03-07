#' Symmetrize a Square Matrix
#'
#' \code{symmetry} transforms a given square matrix into a symmetric matrix.
#'
#' This function symmetrizes the input matrix `a' by averaging it with its
#' transpose, ensuring that the output is symmetric.  Given a square matrix
#' \eqn{A}, the function computes: \deqn{ S = \frac{1}{2} (A + A^T) }
#' This operation ensures that \eqn{S} is symmetric (\eqn{S = S^T}).
#' This is useful in numerical applications where symmetry is required,
#' such as in covariance matrices.
#'
#' @param a A square matrix
#'
#' @return A symmetric square matrix
#'
#' @noRd
#' @examples
#' q <- 3
#' sym.mat <- matrix(rnorm(q * q), nrow = q, ncol = q)
#' sym.mat
#' symmetry(sym.mat)
#'
symmetry <- function(a) {

  # Check if a is a matrix
  if (!is.matrix(a)) {
    stop("a must be a matrix.")
  }

  # Check if a is a square matrix
  if (nrow(a) != ncol(a)) {
    stop("a must be a square matrix.")
  }

  return((a + t(a)) / 2)
}

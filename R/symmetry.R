#' Symmetrize Matrix
#'
#' \code{symmetry} transforms a given square matrix into a symmetric matrix.
#'
#' This function symmetrizes the input matrix `a' by averaging it with its
#' transpose.
#'
#' @param a A square matrix
#'
#' @return A symmetric square matrix
#'
#' @noRd
#' @examples
#' sym.mat <- matrix(rnorm(9), nrow = 3, ncol = 3)
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

  (a + t(a)) / 2
}

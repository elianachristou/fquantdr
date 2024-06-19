#' Symmetrize Matrix
#'
#' \code{symmetry} makes a matrix symmetrical while retaining the most
#' properties possible of the original matrix.
#'
#' The function finds the average between a matrix and its transpose to create
#' a symmetrical matrix while retaining properties such as the mean.
#'
#' @param a The matrix to symmetrize
#'
#' @return A symmetric matrix based on the matrix \code{a}
#'
#' @noRd
#' @examples
#' testMatrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3)
#' symmetry(testMatrix)
#'
symmetry <- function(a) {
  if (!is.matrix(a)) {
    stop("a must be a matrix.")
  }
  if (nrow(a) != ncol(a)) {
    stop("a must be a square matrix.")
  }
  (a + t(a)) / 2
}

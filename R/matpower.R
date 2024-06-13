#' Transformation of a matrix

#' matpower: Power of a matrix

#' \{code} symmetries a matrix using eigen decomposition

#' This function takes in a matrix and an alpha value and proceeds to
#' transform the matrix by computing its eignvectors and eigenvalues, creating a
#' new matrix based on the transpose of a, and then transforming the matrix by
#' creating new matrixs using its eignvectors and eigenvalues and multiply them
#' by themselves.


#' @param a input square matrix
#' @param alpha an exponent to raise the matrix
#'
#' @return returns a matrix that has been raised to a certain power based on its
#' eigenvalues and eigenvectors

#' @noRd

#' @examples
#' matr <- matrix(c(6, 4, 8, 2, 5, 9), nrow = 2, ncol = 3)
#' alpha <- 2
#' matpower(matr, alpha)


matpower <- function(a, alpha) {
  # Computes the transpose of a matrix and then adds it to the original
  # matrix. The new matrix is then divided by two and then rounded to the 5th
  # decimal place.
  a <- round((a + t(a)) / 2, 5)
  # tmp is assigned the eigenvalues and eigenvectors of a
  tmp <- eigen(a)
  # Returns the eigenvectors being multiplied by a diagonal matrix of the
  # eigenvalues where all values on the diagonal have been raised to the power of
  # alpha. This is then multiplied by the transpose of the eigenvectors
  return(tmp$vectors %*% diag((tmp$values)^alpha) %*%
           t(tmp$vectors))}


#' Transformation of a matrix

#' mppower: Moore-Penrose type power

#' \{code} Reconstructs matrix on the basis of significant eigenvalues

#' This code takes an input of a matrix and computes it eigenvalues and
#' eigenvectors. Once computed, it transforms the original matrix by selecting
#' the first m eigen vectors from the matrix eig$vectors and another m eigenvalues
#' from eig$values and then creates a diagonal matrix where the values on the
#' diagonal are the eigenvalues being raised to a power. Matrix multiplication
#' is then preformed to get a transformed matrix that captures the important
#' characteristics of the original matrix.

#' @param matrix an input of a matrix
#' @param power the power we will raise the matrix of eigenvalues to
#' @param ignore a threshold value that denotes to what scale we can ignore the
#' significant vs insignificant eigenvalues
#'
#' @return returns a transformed matrix on the basis of eigenvalues and
#' their corresponding eigenvectors

#' @noRd


#' References

#' @examples
#' matrix1 <- matrix(c(3, 4, 8, 2, 5, 1), nrow = 2, ncol = 3)
#' power <- 3
#' ignore <- .3
#' mmpower(matrix1, power, ignore)



mppower <- function(matrix, power, ignore) {
  # Assigns eig the eigenvectors and eigenvalues of a symmetric matrix
  eig <- eigen(matrix, sym = T)
  # Assignes eval a vector of the eigenvalues
  eval <- eig$values
  # Assigns evec a matrix of the eigenvectors
  evec <- eig$vectors
  # Counts the number of values that have an abs greather than ignore
  m <- length(eval[abs(eval) > ignore])
  # reconstructs the matrix based on significant eigenvalues
  tmp <- evec[, 1:m] %*% diag(eval[1:m]^power) %*% t(evec[, 1:m])
  return(tmp)
}

#' Matrix power transformation
#'
#' \code{matpower} computes the power of a square matrix using eigen
#'     decomposition.
#'
#' This function takes a square matrix \code{a} and an exponent \code{alpha}
#'     and transforms the matrix by performing eigen decomposition.  It
#'     computes the eigenvectors and eigenvalues of the matrix, raises the
#'     eigenvalues to the power of \code{alpha}, and then reconstructs the
#'     matrix using the transformed eigenvectors and adjusted eigenvalues.
#'     The resulting matrix is symmetric (rounded to five decimal places).
#'
#' @param a The input square matrix.
#' @param alpha The exponent to which the matrix is raised.
#'
#' @return The matrix raised to the power of \code{alpha} based on its
#'     eigenvalues and eigenvectors.
#'
#' @noRd
#' @examples
#' mat <- matrix(c(6, 4, 8, 2, 5, 9), nrow = 2, ncol = 3)
#' alpha <- 2
#' matpower(mat, alpha)
#'
matpower <- function(a, alpha) {
  # Checks to make sure matrix returned a matrix
  if (!is.matrix(a)) {
    stop("The first input has not returned a matrix.")
  }
  # Checks if alpha is a positive integer
  if (!is.numeric(alpha) || alpha <= 0 || alpha != as.integer(alpha)) {
    stop("The second value has to be an integer exponenet. ")
  }
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


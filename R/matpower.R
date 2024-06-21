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
  # Checks if 'a' is a matrix
  if (!is.matrix(a)) {
    stop("The input 'a' must be a matrix.")
  }

  # Checks is 'a' is a square matrix
  if (dim(a)[1] != dim(a)[2]) {
    stop("The input 'a' must be a square matrix.")
  }

  # Checks if 'alpha' is a positive integer
  if (!is.numeric(alpha) | alpha <= 0 | alpha != round(alpha)) {
    stop("The exponent 'alpha' must be a positive integer. ")
  }

  # Symmetrize the input matrix 'a'
  a <- round((a + t(a)) / 2, 5)

  # Compute eigen decomposition
  tmp <- eigen(a)

  # Calculate the matrix raised to the power of 'alphaa' based on eigen
  # decomposition
  result <- tmp$vectors %*% diag((tmp$values)^alpha) %*%
           t(tmp$vectors)
  result
}


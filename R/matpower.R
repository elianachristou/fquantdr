#' Matrix power transformation
#'
#' \code{matpower} computes the power of a square matrix using its eigen
#'     decomposition.
#'
#' This function takes a square matrix \code{a} and an exponent \code{alpha}
#' and transforms the matrix by performing eigen decomposition.  It
#' computes the eigenvectors and eigenvalues of the matrix, raises the
#' eigenvalues to the power of \code{alpha}, and then reconstructs the
#' matrix using the transformed eigenvectors and adjusted eigenvalues.
#' The resulting matrix is symmetric (rounded to five decimal places).
#' If \code{alpha} = 1, the function simply returns \code{a}.
#'
#' @param a A square numeric matrix.
#' @param alpha The exponent to which the matrix is raised.
#'
#' @return The square matrix raised to the power of \code{alpha}.
#'
#' @noRd
#' @examples
#' Example 1: Square matrix raised to the power of 2
#' mat <- matrix(c(6, 4, 8, 2, 5, 9, 3, 1, 7), nrow = 3, ncol = 3)
#' alpha <- 2
#' matpower(mat, alpha)
#'
#' Example 2: Computing the square root of a positive definite matrix
#' pos_def_mat <- matrix(c(4, 1, 1, 3), nrow = 2, ncol = 2)
#' matpower(pos_def_mat, 0.5)
#'
matpower <- function(a, alpha) {

  # Check if 'a' is a matrix
  if (!is.matrix(a)) {
    stop("The input 'a' must be a matrix.")
  }

  # Check if 'a' is a square matrix
  if (dim(a)[1] != dim(a)[2]) {
    stop("The input 'a' must be a square matrix.")
  }

  # Check if 'alpha' is a single numeric value
  if (length(alpha) != 1 || !is.numeric(alpha)) {
    stop("The exponent 'alpha' must be a numeric scalar.")
  }

  # Symmetrize the input matrix 'a' for numerical stability
  a <- round((a + t(a)) / 2, 5)

  # Compute eigen decomposition
  tmp <- eigen(a, symmetric = T)

  # Calculate the matrix raised to the power of 'alpha' based on eigen
  # decomposition
  result <- tmp$vectors %*% diag((tmp$values)^alpha) %*%
           t(tmp$vectors)
  return(result)
}


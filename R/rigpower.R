#' Regularized matrix power transformation
#'
#' \code{rigpower} performs a regularized matrix power transformation using
#' eigen decomposition
#'
#' This function takes a square matrix \code{a}, adds a scaled identity
#' matrix based on a regularization parameter \code{rho}, and then raises
#' the resulting matrix to a specified power \code{\alpha}.  This
#' transformation is useful for stabilizing covariace matrices in statistical
#' applications.
#'
#' The regularized transformation is defined as:
#' \deqn{
#' A_{\text{reg}} = A + \rho \lambda_{\max} I
#' }
#' where \eqn{\lambda_{\max}} is the largest eigenvalue of \eqn{A}, and \eqn{I}
#' is the identity matrix. The function then computes:
#' \deqn{
#' A_{\text{reg}}^\alpha
#' }
#'
#' @param a A numeric square matrix.
#' @param alpha A numeric scalar specifying the exponent for the matrix power
#'     transformation.
#' @param rho A numeric scalar specifying the regularization parameter used to
#'     scale.
#'
#' @return A matrix that has been regularized and raised to the specified power.
#'
#' @include matpower.R
#' @noRd
#' @examples
#' # Example: Regularized inverse square root transformation
#' a <- matrix(c(4, 1, 1, 3), nrow = 2, ncol = 2)
#' alpha <- -1 / 2
#' rho <- 0.1
#' a.inv.sqrt <- rigpower(a, alpha, rho)
#' print(a.inv.sqrt)
#'
rigpower <- function(a, alpha, rho) {

  # Check if 'a' is a matrix
  if (!is.matrix(a)) {
    stop("The input 'a' must be a matrix.")
  }

  # Check if 'a' is a square matrix
  if (nrow(a) != ncol(a)) {
    stop("The input 'a' must be a square matrix.")
  }

  # Check if 'alpha' is a numeric scalar
  if (length(alpha) != 1) {
    stop("The exponent 'alpha' must be a single numeric value.")
  }

  # Check if 'alpha' is a numeric scalar
  if (!is.numeric(alpha)) {
    stop("The exponent 'alpha' must be a single numeric value.")
  }

  # Check if 'rho' is a numeric scallar
  if (!is.numeric(rho) || length(rho) != 1) {
    stop("The regularization parameter 'rho' must be a single numeric value.")
  }

  # Eigen decomposition
  eig <- eigen(a, symmetric = T)
  eval <- eig$values

  # Regularize the matrix
  a1 <- a + rho * max(eval) * diag(nrow(a))

  # Apply the matrix power transformation
  tmp <- matpower(a1, alpha)

  return(tmp)
}

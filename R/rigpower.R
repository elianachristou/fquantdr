#' Regularized matrix power transformation
#'
#' \code{rigpower} performs a regularized matrix power transformation using
#' eigen decomposition
#'
#' This function takes a matrix \code{a} and applies a transformation that
#' includes adding a scaled identity matrix to the input matrix and then
#' raising the resulting matrix to a specified power.  The regularization
#' parameter \code{rho} controls the scaling of the identity matrix.
#'
#' @param a The input square matrix.
#' @param alpha The exponent to which the matrix is raised.
#' @param rho A numeric value used as the regularization parameter to scale
#' the identity matrix.
#'
#' @return A matrix that has been regularized and raised to the specified
#' power.
#'
#' @noRd
#' @examples
#' a <- matrix(c(4, 1, 1, 3), nrow = 2, ncol = 2)
#' alpha <- 2
#' rho <- 0.1
#' transformed_matrix <- rigpower(a, alpha, rho)
#' transformed_matrix
#'
rigpower <- function(a, alpha, rho) {

  # Check if 'a' is a matrix
  if (!is.matrix(a)) {
    stop("The input 'a' must be a matrix.")
  }

  # Check is 'a' is a square matrix
  if (dim(a)[1] != dim(a)[2]) {
    stop("The input 'a' must be a square matrix.")
  }

  # Check if 'alpha' is a single number
  if (length(alpha) != 1) {
    stop("The exponent 'alpha' must be a single number.")
  }

  # Check if 'alpha' is numeric
  if (!is.numeric(alpha)) {
    stop("The exponent 'alpha' must be numeric. ")
  }

  # Check if 'rho' is a real number
  if (!is.numeric(rho)) {
    stop("The input 'rho' must be numeric.")
  }

  # define parameters
  p <- nrow(a)

  # Eigen decomposition
  eig <- eigen(a, symmetric = TRUE)
  eval <- Re(eig$values)
  evec <- eig$vectors

  # Regularize the matrix
  a1 <- a + rho * max(eval) * diag(p)

  # Apply the matrix power transformation
  tmp <- matpower(a1, alpha)

  return(tmp)
}

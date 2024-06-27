#' Matrix power transformation
#'
#' \code{mppower} computes the power of a square matrix using eigen
#'     decomposition.
#'
#' This function takes a square matrix \code{a} and an exponent \code{alpha}
#' and transforms the matrix by performing eigen decomposition.  It differs
#' from the function \code{matpower} since it includes a parameter
#' \code{epsilon}, which, if provided, it is added to the diagonal elements
#' to stabilize computations.  Moreover, it includes a parameter \code{ignore},
#' which provides a threshold below which eigenvalues are ignored.
#'
#' @param a The input square matrix.
#' @param alpha The exponent to which the matrix is raised.
#' @param epsilon A nonnegative numeric value added to the diagonal elements
#'     of the matrix to stabilize computations (default is 0).
#' @param ignore A numeric threshold below which eigenvalues are ignored
#'     (default is 10^(-15)).
#'
#' @return The matrix raised to the power of \code{alpha}.

#' @noRd
#' @examples
#' mat <- matrix(c(4, 1, 1, 3), 2, 2)
#' result <- mppower(mat, 2)
#'
mppower <- function(matrix, power, epsilon = 0, ignore = 10^(-15)) {
  # Checks if the matrix input is a matrix
  if (!is.matrix(matrix)){
    stop("The first input must be a matrix.")
  }
  # Checks if power is a positive integer
  #if (!is.numeric(power) || power <= 0 || power != as.integer(power)) {
  #  stop("The second value has to be an integer exponenet. ")
  #}
  # Compatability check for ignore
  if (!is.numeric(ignore)) {
    stop("Must be a real number.")
  }
  # Assigns eig the eigenvectors and eigenvalues of a symmetric matrix
  eig <- eigen(matrix, symmetric = T)
  # Assignes eval a vector of the eigenvalues
  eval <- eig$values
  # Assigns evec a matrix of the eigenvectors
  evec <- eig$vectors
  # Counts the number of values that have an abs greather than ignore
  m <- length(eval[abs(eval) > ignore])
  # reconstructs the matrix based on significant eigenvalues
  tmp <- evec[, 1:m] %*% diag(eval[1:m]^power) %*% t(evec[, 1:m])
  # Checks if the computation produced a matrix
  if (!is.matrix(tmp)) {
    stop("tmp must be a matrix.")
  }
  return(tmp)
}

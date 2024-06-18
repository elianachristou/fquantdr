#' Preforms a power operation on a transformed input matrix

#' \code{rigpower} This function takes in three values: a matrix, a power, and the value
#' rho. Using those values the function computes the eigenvalues and eigenvectors
#' of the matrix. It uses that information to then transform the original matrix
#' using the matpower function.


#' @param matrix a matrix input
#' @param power an exponement that raises the modified matrix
#' @param rho a scaling value that is multiplied by the maximum eigenvector
#'
#' @return returns a modified version of the input matrix that has been altered
#' by an scaled identity matrix
#' @noRd

#' @examples
#' matx <- matrix(c(7, 3, 5, 2, 5, 1), nrow = 2, ncol = 3)
#' power<- 2
#' rho <- .1
#' rigpower(matx, power, rho)

rigpower <- function(matrix, power, rho) {
  # Checks to make sure matrix returned a matrix
  if (!is.matrix(matrix)) {
    stop("The first input has not returned a matrix.")
  }
  # Checks if power is a positive integer
  if (!is.numeric(power) || power <= 0 || power != as.integer(power)) {
    stop("The second value has to be an integer exponenet. ")
  }
  if (!is.numeric(rho)) {
    stop("Must be a real number.")
  }
  # eig is assigned the eigenvalues and vectors of the matrix
  eig <- eigen(matrix)
  # eval is assigned a vector of the eigenvalues
  eval <- eig$values
  # Assigned a matrix of the eigenvectors
  evec <- eig$vectors
  # p is assigned the number of rows in the matrix
  p <- nrow(matrix)
  # Assigns matrix1 the value of the original matrix being incremented by a
  # diagonal matrix that stores the values of the maximum eigenvalue being
  # multiplied by rho
  matrix1 <- matrix + rho * max(eval) * diag(p)
  # Checks if the computation produced a matrix
  if (!is.matrix(matrix1)) {
    stop("matrix1 must be a matrix.")
  }
  # tmp is assigned the value of a transformed matrix from matpower
  tmp <- matpower(matrix1, power)
  # Checks if the computation produced a matrix
  if (!is.matrix(tmp)) {
    stop("tmp must be a matrix.")
  }
  # returns the matrix
  return(tmp) }

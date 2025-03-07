#' Projection Matrix onto the Orthogonal Complement of the Subspace spanned by 1
#'
#' \code{qmat} computes the projection matrix onto the orthogonal complement
#' of the subspace spanned by the vector of ones.
#'
#' This function generates a matrix that removes the mean by projecting
#' onto the orthogonal complement of the subspace spanned by a vector of ones.
#' This is useful in various statistical and mathematical applications where
#' centering or removing the mean is required.  The projection matrix is given by:
#' \deqn{Q = I_n - (1/n) J_n}
#' where:
#' \itemize{
#'   \item \eqn{I_n} is the \eqn{n \times n} identity matrix.
#'   \item \eqn{J_n} is an \eqn{n \times n} matrix of ones.
#' }
#' This matrix is symmetric and idempotent, meaning \eqn{Q^2 = Q}.
#'
#' @param n A positive integer representing the number of rows and columns
#' of the resulting square matrix.
#'
#' @return An \code{n x n} matrix representing the projection onto the
#' orthogonal complement of the subspace spanned by the vector of ones.
#'
#' @noRd
#' @examples
#' # Generate the 3 x 3 projection matrix
#' qmat(3)
#'
qmat <- function(n) {

  # Check if n is a positive integer
  if (n <= 0 | n != floor(n)) {
    stop("n must be a positive integer.")
  }

  # Create the n x n identity matrix and the matrix of ones
  ii <- diag(n)
  jj <- rep(1, n) %*% t(rep(1, n))

  # Return the projection matrix
  return(ii - jj / n)
  }


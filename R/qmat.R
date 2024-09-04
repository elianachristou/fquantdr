#' Projection onto the orthogonal complement of the subspace spanned by 1
#'
#' \code{qmat} computes the projection matrix onto the orthogonal complement
#' of the subspace spanned by the vector of ones.
#'
#' This function generates a matrix that projects onto the orthogonal
#' complement of the subspace spanned by a vector of ones. This is useful
#' in various statistical and mathematical applications where centering or
#' removing the mean is required.
#'
#' @param n A positive integer representing the number of rows and columns
#' of the resulting square matrix.
#'
#' @return A \code{n x n} matrix representing the projection onto the
#' orthogonal complement of the subspace spanned by the vector of ones.
#'
#' @noRd
#' @examples
#' # Generate the 3 x 3 projection matrix
#' proj.matrix3 <- qmat(3)
#' proj.matrix3
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


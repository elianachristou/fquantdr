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
#'
#' @seealso \code{\link{matpower}}
#'
#' @noRd
#' @examples
#' mat <- matrix(c(6, 4, 8, 2, 5, 9, 3, 1, 7), nrow = 3, ncol = 3)
#' alpha <- 2
#' result <- mppower(mat, alpha)
#'
mppower <- function(a, alpha, epsilon = 0, ignore = 10^(-15)) {

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

  # Check if 'epsilon is >= 0'
  if (epsilon < 0){
    stop("epsilon should be a nonnegative real number.")
  }
  # Checks if epsilon is complex
  if (is.complex(epsilon)) {
    stop("epsilon should be a nonnegative real number.")
  }

  # to ensure that the matrix is symmetric
  B <- symmetry(a)

  if (epsilon > 0) {
    iden <- diag(nrow(a))
    B <- B + epsilon * iden
  }

  # eigen decomposition
  eig <- eigen(B, symmetric = T)
  eval <- eig$values
  evec <- eig$vectors

  if(!is.numeric(ignore)) {
    stop("the third input, ignore, must be an numeric input")
  }

  # Find indices where eigenvalues are greater than ignore threshold
  m <- length(eval[eval > ignore])
  tmp <- evec[, 1:m] %*% diag(eval[1:m]^alpha) %*% t(evec[, 1:m])
  return(tmp)
}

#' Computes a transformed identity matrix

#' qmat: Q = I - J/n matrix

#' \{code} Creates a variation of the identity matrix such that each row sums
#' zero

#' This function takes a value for n and uses it to create an identity
#' matrix with n columns and n rows. It then creates another matrix with a size
#' of n x n where every entry is assigned the value 1. The function then divides
#' the all-ones matrix by n and then subtracts the identity matrix from the new
#' matrix.

#' @param n a value being assigned at the input area

#' @return gives back a I - J/n matrix

#' @noRd

#' @examples
#' #number of columns
#' n <- 3
#' # An 3 x 3 identity matrix
#' ii = diag(3)
#' # A 3 x 3 all-ones matrix
#' jj <- rep(1, 3) %*% t(rep(1, 3))
#' # returns a new matrix
#' qmat(3)

qmat <- function(n) {
  # Creates a n x n identity matrix
  ii <- diag(n)
  # Creates a n x n matrix where every value is assigned 1
  jj <- rep(1, n) %*% t(rep(1, n))
  # divides the all-ones matrix by n and subtracts the new matrix from the identiy matrix
  return(ii - jj / n) }


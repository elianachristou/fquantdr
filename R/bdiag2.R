#' Block diagonal matrix
#'
#' \code{bdiag2} constructs a block diagonal matrix from two input matrices.
#'
#' This function takes two matrices and arranges them in block diagonal form,
#' meaning `M` occupies the top-left block and `N` occupies the bottom-right block,
#' with zeros filling the remaining space.  The resulting matrix has dimensions
#' `(nrow(M) + nrow(N)) \times (ncol(M) + ncol(N))`.
#'
#' @param M,N Two numeric matrices
#'
#' @return A block diagonal matrix combining `M` and `N`
#'
#' @noRd
#' @examples
#' M <- matrix(1:4, 2, 2)
#' N <- matrix(5:8, 2, 2)
#' bdiag2(M, N)
#'
bdiag2 <- function(M, N) {
  # Get dimensions of input matries
  nc <- ncol(M)
  nr <- nrow(M)
  nc2 <- ncol(N)
  nr2 <- nrow(N)

  # Initialize an empty matrix with appropriate size
  out <- matrix(0, nr + nr2, nc + nc2)
  # Pllace M in the top-left block
  out[(1:nr), (1:nc)] <- M
  # Place N in the bottom-right block
  out[((nr + 1):(nr + nr2)), ((nc + 1):(nc + nc2))] <- N
  return(out)
}

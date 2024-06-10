#' Symmetrize Matrix
#'
#' \code{symmetry} makes a matrix symmetrical while retaining the most
#' properties possible of the original matrix.
#'
#' The function finds the average between a matrix and its transpose to create
#' a symmetrical matrix while retaining properties such as the mean and
#'
#' @param a
#'
#' @return
#' @export
#'
#' @examples
symmetry <- function(a) {
  (a + t(a)) / 2
}

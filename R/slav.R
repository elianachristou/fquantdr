#' Group-wise mean calculation
#'
#' \code{slav} computes the mean values of \code{x} (row-wise) within each
#' group determined by the \code{H} slices.
#'
#' This function calculates the mean of \code{x} (row-wise) for each group
#' defined by the unique elements of the sliced response \code{y}.  The
#' result is a matrix containing the mean values for each group.  This is
#' essential for applying sliced inverse regression (SIR) of Li (1991) and
#' its functional counterpart of Ferr\'e and Yao (2003).
#'
#' @param x A \code{n x p} numeric matrix, where rows represent
#'     observations and columns represent variables.  If the user is using
#'     functional predictors, then \code{x} consists of the coefficients
#'     of the functional object.
#' @param y A numeric vector.
#' @param H The number of slices.
#'
#' @return A \code{n x p} matrix, where rows represent the mean values of
#'     \code{x} as determined by slicing the response \code{y}.
#'
#' @references Li, K.-C. (1991) Sliced Inverse Regression for Dimension Reduction.
#' \emph{Journal of the American Statistical Association}, 86(414), 316-327.
#'
#' Ferr\'e, L, and Yao, F. (2003) Function Sliced Inverse Regression
#' Analysis. \emph{Statistics}, 37(6), 475-488.
#'
#' @seealso \code{\link{discretize}}
#'
#' @noRd
#' @examples
#' x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
#'   nrow = 4, byrow= TRUE)
#' y <- c(1, 2, 3, 2)
#' slav(x, y, 3)
#'
slav <- function(x, y, H) {

  # Check to ensure x is a matrix
  if (!is.matrix(x)){
    stop("The input 'x' must be a matrix")
  }

  if (!is.vector(y)) {
    stop("The input 'y' must be a vector.")
  }

  if (length(y) != dim(x)[1]) {
    stop("The number of observations of y must be the same as the number
         of rows of x.")
  }

  if (!is.numeric(H) | H <= 0 | H != floor(H)) {
    stop("H must be a positive integer.")
  }

  for (i in y) {
    if(!(i %in% 1:H)) {
      stop(paste("y must be discretized such that all values in y are in",
                 "yunit. The value", i, "is not in yunit."))
    }
  }

  # Define the parameters
  n <- nrow(x)
  p <- ncol(x)
  yunit <- 1:H

  # Initialize matrix to store group means
  xgy <- matrix(0, H, p)

  # Compute the mean of x for each group defined by yunit
  for(i in yunit) {
    predictors <- x[y == yunit[i], ]
    # If vector, convert to horizontal matrix
    if (is.vector(predictors)) {
      predictors <- t(as.matrix(predictors))
    }
    xgy[i, ] <- apply(predictors, 2, mean)
  }

  xgy
}

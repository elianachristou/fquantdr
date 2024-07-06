#' Group-wise mean calculation
#'
#' \code{slav} computes the mean of the rows of \code{x} within each group
#' determined by the \code{H} slices.
#'
#' This function calculates the mean of the rows of \code{x} for each group
#' defined by the unique elements of the sliced response \code{y}.  The
#' resulting matrix contains the mean values for each group.  This is
#' essential for applying sliced inverse regression (SIR) of Li (1991) and
#' its functional counterpart of Ferr\'e and Yao (2003).
#'
#' @param x A \code{n x p} numeric matrix, where rows represent
#'     observations and columns represent variables.  If the user is using
#'     functional predictors, then \code{x} consists of the coefficients
#'     of the functional object and is of \code{n x (p * nbasis)} dimension.
#' @param y A numeric vector.
#' @param H The number of slices.
#'
#' @return A \code{H x p} matrix, where rows represent the mean values of
#'     \code{x} as determined by slicing the response \code{y}.
#'
#' @references Li, K.-C. (1991) Sliced Inverse Regression for Dimension
#' Reduction. \emph{Journal of the American Statistical Association}, 86(414),
#' 316-327.
#'
#' Ferr\'e, L, and Yao, F. (2003) Function Sliced Inverse Regression Analysis.
#' \emph{Statistics}, 37(6), 475-488.
#'
#' @seealso \code{\link{discretize}}, \code{\link{slprob}}
#'
#' @noRd
#' @examples
#' n <- 20
#' p <- 5
#' H <- 4
#' x <- matrix(rnorm(100), nrow = n, ncol = p)
#' y <- rnorm(n)
#' slav(x, y, H)
#'
slav <- function(x, y, H) {

  # Check if x is a matrix
  if (!is.matrix(x)){
    stop("The input 'x' must be a matrix.")
  }

  # Check if y is a vector
  if (!is.vector(y)) {
    stop("The input 'y' must be a vector.")
  }

  # Check if x and y have the same number of observations
  if (length(y) != dim(x)[1]) {
    stop("The number of observations of y must be the same as the number
         of rows of x.")
  }

  # Check if H is a positive integer
  if (H <= 0 | H != floor(H)) {
    stop("H must be a positive integer.")
  }

  # Define the parameters
  n <- nrow(x)
  p <- ncol(x)
  yunit <- 1:H
  ydis <- discretize (y, H)

  # Initialize matrix to store group means
  xgy <- matrix(0, H, p)

  # Compute the mean of x for each group defined by yunit
  for(i in yunit) {
    predictors <- x[ydis == yunit[i], ]
    # If vector, convert to horizontal matrix
    if (is.vector(predictors)) {
      predictors <- t(as.matrix(predictors))
    }
    xgy[i, ] <- apply(predictors, 2, mean)
  }

  xgy
}

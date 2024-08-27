#' Group-wise mean calculation
#'
#' \code{slav} computes the row-wise mean of \code{x} within each group
#' determined by the \code{H} slices.
#'
#' This function calculates the row-wise mean of \code{x} for each group
#' defined by the unique elements of the sliced response \code{y}.  The
#' resulting matrix contains the mean values for each group.  This is
#' essential for applying sliced inverse regression (SIR) of Li (1991) and
#' its functional counterpart of Ferr\'e and Yao (2003).
#'
#' @param x A \code{n x p} numeric matrix, where rows represent
#'     observations and columns represent variables.  If the user is using
#'     functional predictors, then \code{x} consists of the coefficients
#'     of the functional object and is of \code{n x (p * nbasis)} dimension.
#' @param ydis A discrete numeric vector that contains labels that define the
#' slice each value of a vector `y` is in. This can be obtained using the
#' `discretize` function.
#' @param H The number of slices.
#'
#' @return A \code{H x p} matrix, where rows represent the mean values of
#'     \code{x} as determined by the discretize response \code{ydis}. If
#'     \code{x} consists of the coefficients of a functional object, then
#'     the output is a \code{H x (p * nbasis)} matrix.
#'
#' @references Li, K.-C. (1991) Sliced Inverse Regression for Dimension
#' Reduction. \emph{Journal of the American Statistical Association}, 86(414),
#' 316-327.
#'
#' Ferr&#233;, L, and Yao, F. (2003) Function Sliced Inverse Regression Analysis.
#' \emph{Statistics}, 37(6), 475-488.
#'
#' @seealso \code{\link{discretize}}, \code{\link{slprob}}
#'
#' @noRd
#' @examples
#' # Example 1
#' # set the parameters
#' n <- 100
#' p <- 5
#' H <- 10
#' # define the variables
#' x <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' beta <- c(1, 1, 0, 0, 0)
#' y <- c(x %*% beta + rnorm(n))
#' ydis <- discretize(y, H)
#' slav(x, ydis, H)
#'
#' # Example 2
#' # set the parameters
#' n <- 100
#' p <- 5
#' nt <- 101
#' nbasis <- 4
#' H <- 10
#' time <- seq(0, 1, length.out = nt)
#' eta <- matrix(stats::rnorm(n * p * nbasis), nrow = n,
#'     ncol = p * nbasis)
#' # Generate the functional data
#' result <- fundata(n, p, nbasis, time, eta)
#' Xc <- result$xc
#' P <- eigen(stats::cov(eta))$vectors
#' mfpca.scores <- eta %*% P
#' # Generate the model
#' error <- rnorm(n)
#' y <- 3 * mfpca.scores[, 1] + error
#' ydis <- discretize(y, H)
#' # define the coefficients for the functional predictors
#' databasis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
#' xfd.coef <- numeric()
#' for (k in 1:p) {
#'  xfdk <- fda::smooth.basis(seq(0, 1, length.out = nt),
#'                            t(Xc[, , k]), databasis)$fd
#'  xfdk <- fda::center.fd(xfdk)
#'  xk.coef <- t(xfdk$coef)
#'  xfd.coef <- cbind(xfd.coef, xk.coef)
#' }
#' # Calculate the row-wise means
#' slav(xfd.coef, ydis, H)
#'
slav <- function(x, ydis, H) {

  # Check if x is a matrix
  if (!is.matrix(x)){
    stop("The input 'x' must be a matrix.")
  }

  # Check if ydis is a vector
  if (!is.vector(ydis)) {
    stop("The input 'ydis' must be a vector.")
  }

  # Check if x and ydis have the same number of observations
  if (length(ydis) != dim(x)[1]) {
    stop("The number of observations of ydis must be the same as the number
         of rows of x.")
  }

  # Check if H is a positive integer
  if (H <= 0 | H != floor(H)) {
    stop("H must be a positive integer.")
  }

  # Check if H is one number
  if (length(H) > 1) {
    stop("H must be one number.")
  }

  # Check if H is greater than 1 and less than the length of ydis
  if (H == 1 | H >= length(ydis)) {
    stop("H must be an integer that is at least 2 and less than the
         length of ydis.")
  }

  # Define the parameters
  n <- nrow(x)
  p <- ncol(x)
  yunit <- 1:H

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

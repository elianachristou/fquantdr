#' Functional Principal Component Analysis (FPCA)
#'
#' \code{fpca} computes FPCA by extracting the coefficient matrix from the
#' functional data object, constructing the appropriate penalty matrix based
#' on the specified basis type, and applying eigen decomposition to obtain
#' functional principal components.
#'
#' This function performs FPCA on a given function object and returns the
#' eigenvalues, transformation matrix, and principal component scores of the
#' functional covariance matrix.  The method supports both univariate
#' (\code{p = 1}) and multivariaate (\code{p > 1}) functional predictors.
#'
#' @param ftn A functional data object containing coefficients and basis
#'     information.  The coefficients is a 3D array (nbasis x n x p), where `nbasis`
#'     is the number of basis, `n` is the sample size, and `p` is the number of
#'     functional variables.  For the basis, options are 'bspline' or 'fourier'.
#'
#' @return \code{fpca} returns a list containing:
#'    \itemize{
#'        \item \code{pred}: Principal component scores.
#'        \item \code{eval}: Eigenvalues of the covariance operator.
#'        \item \code{mat}: Transformation matrix for principal components.
#'    }
#'
#' @examples
#' # set the parameters
#' n <- 100
#' p <- 5
#' nbasis <- 4
#' nt <- 100
#' tt <- seq(0, 1, length.out = nt)
#' eta.mat <- mvtnorm::rmvnorm(n, mean = rep(0, p * nbasis))
#' eta <- array(eta.mat, dim = c(nbasis, n, p))
#' # Generate the functional data
#' data <- fundata(n, p, nbasis, tt, 'bspline', eta)
#' pca.out <- fpca(list(coef = data$xcoefs, basis = data$basis), 'bspline')
#' pca.out$eval
#' mfpca.scores <- pca.out$pred
#'
#' @export
fpca <- function(ftn) {

  # Check if ftn is a valid functional object
  if (!("coef" %in% names(ftn)) || !("basis" %in% names(ftn))) {
    stop("'ftn' must be a functional object with 'coef' and 'basis' components.")
  }

  # Check if basisname is valid
  if (!ftn$basis$type %in% c("bspline", "fourier")) {
    stop("'basisname' must be either 'bspline' or 'fourier'.")
  }

  # Check if coef array has expected dimensions
  if (!is.array(ftn$coef) || length(dim(ftn$coef)) < 2) {
    stop("'ftn$coef' must be a multidimensional array, where the first dimension
         represent the number of basis functions, the second dimension represent
         the number of observations, and the third dimension represent the number
         of variables.  If the third dimension is not specified, then it will be
         set to 1, treating it as a univariate case.")
  }

  # Extract coefficient array from input 'ftn' and determine its dimensions
  temp <- ftn$coef
  n <- dim(temp)[2]
  p <- dim(temp)[3]
  nbasis <- dim(temp)[1]

  # If 'p' is NA, set it to 1 (for univariate case)
  if(is.na(p)) p <- 1
  xcoef <- temp

  # Extract basis information from input
  basis <- ftn$basis
  basisname <- basis$type

  # Generate basis penalty matrix based on specified basis type
  if(basisname == 'bspline') {
    GB <- fda::bsplinepen(basis, Lfdobj = 0)
  }
  else if(basisname == 'fourier') {
    GB <- fda::fourierpen(basis, Lfdobj = 0)
  }

  # Create centering matrix for FPCA
  one <- matrix(1, n, 1)
  Q <- diag(n) - one %*% t(one) / n

  # Compute the square root of the penalty matrix
  B.half <- matpower(GB, 0.5)

  # Univariate case
  if(p==1) {
    # Compute covariance matrix in transformed space
    Sigma <- B.half %*% xcoef %*% Q %*% t(xcoef) %*% B.half / n
    # Perform eigen decomposition of covariacne matrix
    egn <- eigen(Sigma, symmetric = TRUE)
    B.inv.half <- matpower(GB, -0.5)
    # Compute predicited principal component scores
    pred <- Q %*% t(xcoef) %*% GB %*% B.inv.half %*% egn$vec
    out <- list(pred = pred, eval = egn$val, mat = B.inv.half %*% egn$vec)
  } else if(p > 1) {
    # Multivariate case
    # Compute transformed coefficient matrix
    M.half <- B.half %*% xcoef[, , 1] %*% Q
    B.inv.half <- matpower(GB, -0.5)
    D.inv.half <- B.inv.half
    # Compute transformed matrix for all variables
    BB <- GB %*% xcoef[, , 1] %*% Q
    for(j in 2:p) {
      D.inv.half <- as.matrix(bdiag2(D.inv.half, B.inv.half))
      M.half <- rbind(M.half, B.half %*% xcoef[, , j] %*% Q)
      BB <- rbind(BB, GB %*% xcoef[, , j] %*% Q)
    }
    # Compute covariance matrix in transformed space
    M.half <- M.half %*% t(M.half) / n
    egn <- eigen(M.half, symmetric = TRUE)
    # Compute predicted principall component scores
    pred <- t(BB) %*% D.inv.half %*% egn$vec
    out <- list(pred = pred, eval = egn$val, mat = D.inv.half %*% egn$vec)
  }
  return(out)
}

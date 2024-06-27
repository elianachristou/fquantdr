#' Functional Sliced Inverse Regression
#'
#' \code{mfsir} performs dimension reduction for functional data and
#' provides the new estimated predictors.
#'
#' This function performs functional sliced inverse regression (FSIR),
#' introduced by Ferr\'e and Yao (2003). Specifically, the authors proved
#' that \eqn{E(X|Y) - E(X)} belongs to \eqn{\Sigma_{XX} S_{Y|X}}, where
#' \eqn{S_{Y|X}} denotes the functional central subspace.
#'
#' @param X A \code{n x t x p} array, where n is the number of observations,
#'     t is the number of time points, and p is the number of predictor
#'     variables.
#' @param y A numeric vector of length \code{n} representing the response
#'     variable.
#' @param H The number of slices.
#' @param nbasis The number of basis functions for smoothing functional data.
#'
#' @return \code{mfsir} computes the new sufficient predictors and returns
#' \item{xcoef}{The coefficients of the smoothed functional data.}
#' \item{eigvalues}{The eigenvalues resulting from the eigenvalue decomposition
#'    of the matrix that gives the sufficient predictors.}
#' \item{phi}{The eigenvectors resulting from the eigenvalue decomposition
#'     of the matrix that gives the sufficient predictors.}
#' \item{gx}{The block diagonal Gram matrix.}
#' \item{betas}{The estimated coefficients that form the sufficient
#'     predictors.}
#' \item{sufpred}{The estimated sufficient predictors of the functional
#'     central subspace.}
#'
#' @examples
#' # set the parameters
#' n <- 100
#' p <- 5
#' nt <- 101
#' nbasis <- 4
#' H <- 10
#' time <- seq(0, 1, length.out = nt)
#' eta <- matrix(stats::rnorm(n * p * nbasis), nrow = n,
#'     ncol = p * nbasis)
#'
#' # Generate the functional data
#' result <- fundata(n, p, nbasis, time, eta)
#' X <- result$cg
#' P <- eigen(stats::cov(eta))$vectors
#' mfpca.scores <- eta %*% P
#'
#' # Generate the model
#' error <- rnorm(n)
#' y <- 3 * mfpca.scores[, 1] + error
#'
#' # Run mfsir
#' result <- mfsir(X, y, H, nbasis)
#' result$sufpred
#' # Plot the first sufficient predictor against the true one
#' plot(result$sufpred[, 1], mfpca.scores[, 1], xlab = 'First
#'     Sufficient Predictor', ylab = 'True Predictor')
#'
#' @noRd
mfsir <- function(X, y, H, nbasis) {

  # Check if X is a 3-dimensional array
  if (length(dim(X)) != 3) {
    stop("X must be a 3-dimensional array, where the first dimension
         represent the number of observations, the second dimension
         represent the number of time points, and the third dimension
         reepresent the number of predictor variables.")
  }

  # Check if y is a vector of appropriate length
  if (!is.vector(y) || length(y) != dim(Xc)[1]) {
    stop("y must be a vector of length equal to n, the number of
    observations in X")
  }

  # Check if H is a positive integer
  if (!is.numeric(H) || H <= 0 || H != floor(H)) {
    stop("H must be a positive integer")
  }

  # Check if nbasis is a positive integer
  if (!is.numeric(nbasis) || nbasis <= 0 || nbasis != floor(nbasis)) {
    stop("nbasis must be a positive integer")
  }

  # define the parameters
  n <- dim(X)[1]
  nt <- dim(X)[2]
  p <- dim(X)[3]

  # center X
  Xc <- array(0, dim = c(n, nt, p))
  for(k in 1:p) {
  Xc[, , k] <- X[, , k] - matrix(rep(apply(X[, , k], 2, mean), n),
                                 nrow = n, byrow = T)
  }

  # Create a B-spline basis for smoothing
  databasis <- fda::create.bspline.basis(rangeval = c(0, 1),
                                         nbasis = nbasis)

  # Calculate the coefficients
  xfd.coef <- numeric()
  for (k in 1:p) {
    xfdk <- fda::smooth.basis(seq(0, 1, length.out = nt),
                              t(Xc[, , k]), databasis)$fd
    xfdk <- fda::center.fd(xfdk)
    xk.coef <- t(xfdk$coef)
    xfd.coef <- cbind(xfd.coef, xk.coef)
  }

  #compute the block diagonal gram matrix gx
  gi <- gramatrix(nbasis, databasis)
  gx <- matrix(0, nrow = p * (nbasis), ncol = p * (nbasis))
  for (i in 1:p) {
    index <- ((i - 1) * (nbasis) + 1):(i * (nbasis))
    gx[index, index] <- gi
  }

  # Compute the covariance matrix
  sigmaxx <- (t(xfd.coef) %*% xfd.coef) / n

  # Compute the square root of the Gram matrix
  gx.sqrt <- mppower(gx, 1/2, 10e-7)

  # Compute the target matrix and its inverse square root
  store <- gx.sqrt %*% sigmaxx %*% gx.sqrt
  store.inv.sqrt <- rigpower(store, -1/2, 0.05)

  # Compute the mean of the coefficients
  mucoef <- apply(xfd.coef, 2, mean)

  # Discretize the response variable into H slices and compute the slice
  # probabilities and the means within each slice
  ydis <- discretize (y, H)
  prob <- slprob(ydis, H)
  avg.slav <- slav(xfd.coef, ydis, H)

  # Calculate the target matrix for the eigenvalue decomposition
  Lambda1 <- matrix(0, p * nbasis, p * nbasis)
  for(i in 1:H) {
    Lambda1 <- Lambda1 + prob[i] * avg.slav[i, ] %*% t(avg.slav[i, ])
    }

  # Compute the cross-covariance matrix
  sigmaxxy <- Lambda1 - mucoef %*% (t(mucoef))

  # Compute the transformation matrix
  mat <- store.inv.sqrt %*% gx.sqrt

  # Perform eigenvalue decomposition
  out <- mat %*% sigmaxxy %*% t(mat)
  out <- symmetry(out)
  eigendecom <- eigen(out)
  eigvalues <- eigendecom$values
  phi <- eigendecom$vec

  # Compute the inverse square root of the Gram matrix
  gx.inv.sqrt <- rigpower(gx, -1/2, 0.05)

  # Compute the regression coefficients and the sufficient predictors
  betas <- gx.inv.sqrt %*% store.inv.sqrt %*% phi
  vv <- xfd.coef %*% gx %*% betas

  return(list(xcoef = xfd.coef, eigvalues = eigvalues, phi = phi,
              gx = gx, betas = betas, sufpred = vv))
}

#' Functional Sliced Inverse Regression
#'
#' \code{mfsir} performs dimension reduction for functional data and
#' provides the new estimated predictors.
#'
#' This function implements functional sliced inverse regression (FSIR),
#' introduced by Ferr&#233; and Yao (2003), for scalar-on-function
#' regression problems.  The method seeks to identify a low-dimensional
#' subspace that captures the relationship between a scalar response and
#' functional predictors.  The authors proved that \eqn{E(X|Y) - E(X)}
#' belongs to \eqn{\Sigma_{XX} S_{Y|X}}, where \eqn{S_{Y|X}} denotes the
#' functional central subspace.
#'
#' For \eqn{i=1, \dots, p}, let \eqn{\mathcal{H}_i} be a separable
#' Hilbert space of real-valued functions on \eqn{T}, a bounded closed interval
#' in \eqn{\mathbb{R}}.  Let \eqn{Y: \Omega \rightarrow \mathbb{R}} a univariate
#' response and \eqn{X=(X^1, \dots, X^p): \Omega \rightarrow
#' \bigoplus_{i=1}^p \mathcal{H}_i} a random element.  Dimension reduction
#' techniques aim at finding functions \eqn{\beta_1, \dots, \beta_d} in
#' \eqn{\bigoplus_{i=1}^p \mathcal{H}_i}, such that
#' \deqn{Y = g(\langle \beta_1, X \rangle_{\bigoplus \mathcal{H}}, \dots,
#' \langle \beta_d, X \rangle_{\bigoplus \mathcal{H}}, \epsilon),} where
#' \eqn{g} is an arbitrary unknown function on \eqn{\mathbb{R}^{d+1}}, and
#' \eqn{\epsilon} is independent of \eqn{X}.  This implies that \eqn{Y} and
#' \eqn{X} are independent given \eqn{\langle \beta_1,
#' X \rangle_{\bigoplus \mathcal{H}}, \dots,
#' \langle \beta_d, X \rangle_{\bigoplus \mathcal{H}}} and that the \eqn{p}-
#' dimensional predictor \eqn{X} can be replaced with the \eqn{d}-dimensional
#' predictor \eqn{\langle \beta_1, X \rangle_{\bigoplus \mathcal{H}}, \dots,
#' \langle \beta_d, X \rangle_{\bigoplus \mathcal{H}}}.
#'
#' The functions \eqn{\beta_1, \dots, \beta_d} are called the \emph{functional
#' dimension reduction directions} and the subspace spanned by \eqn{\beta_1,
#' \dots, \beta_d} is called the \emph{functional dimension reduction
#' subspace}.  The smallest functional dimension reduction subspace is called
#' the \emph{functional central subspace} and is denoted by \eqn{S_{Y|X}}.
#'
#' @param X A 3-dimensional array (\code{n x nt x p}), where n is the number
#'     of observations, nt is the number of time points, and p is the number
#'     of predictor variables.
#' @param y A numeric vector of length \code{n} representing the scalar response
#'     variable.
#' @param H An integer specifying the number of slices for the response variable;
#'     must be at least 2 and less than \code{n}.
#' @param nbasis An integer specifying the number of basis functions for functional
#'     smoothing.  Currently, the only option is to use B-spline basis. Default
#'     value is 4.
#' @param norder An integer specifying the order of B-splines, which is one
#'     higher than their degree.  The default of 4 gives cubic splines.
#'
#' @return \code{mfsir} computes the new sufficient predictors and returns
#' \item{xcoef}{A \code{n x (p * nbasis)} matrix of smoothed and centered
#' coefficients for the functional predictors.}
#' \item{eigvalues}{A vector of eigenvalues obtained from the eigenvalue
#' decomposition of the matrix of interest that is calculated during the
#' dimension reduction process.}
#' \item{phi}{A \code{(p * nbasis) x (p * nbasis)} matrix of eigenvectors
#'      resulting from the eigenvalue decomposition of the matrix of interest
#'      that is calculated during the dimension reduction process.}
#' \item{gx}{The \code{(p * nbasis) x (p * nbasis)} block diagonal Gram
#'     matrix of the B-spline basis functions.}
#' \item{betas}{The \code{(p * nbasis) x (p * nbasis)} matrix of the
#'     coordinates of \eqn{\beta_1, \dots, \beta_d} resulting from the
#'     coordinate representation on the B-spline basis functions.}
#' \item{sufpred}{A \code{n x (p * nbasis)} matrix of estimated sufficient
#' predictors of the functional central subspace.}
#'
#' @references Ferr&#233;, L, and Yao, F. (2003) Function Sliced Inverse Regression
#' Analysis. \emph{Statistics}, 37(6), 475-488.
#'
#' @examples
#' # set the parameters
#' n <- 100
#' p <- 5
#' nt <- 100
#' nbasis <- 4
#' H <- 10
#' tt <- seq(0, 1, length.out = nt)
#' eta.mat <- mvtnorm::rmvnorm(n, mean = rep(0, p * nbasis))
#' eta <- array(eta.mat, dim = c(nbasis, n, p))
#' # Generate the functional data
#' data <- fundata(n, p, nbasis, tt, 'fourier', eta)
#' Xc <- data$Xc
#' # Generate the model
#' error <- rnorm(n)
#' y <- 3 * data$mfpca.scores[, 1] + error
#' # Run mfsir
#' result <- mfsir(Xc, y, H, nbasis)
#' result$sufpred
#' # Plot the first sufficient predictor against the true one
#' plot(result$sufpred[, 1], data$mfpca.scores[, 1], xlab = 'First
#'     Sufficient Predictor', ylab = 'First true Predictor', pch = 20)
#' # Calculate the correlation between the estimated and true predictors
#' mcorr(result$sufpred[, 1], data$mfpca.scores[, 1])
#'
#' @export
mfsir <- function(X, y, H, nbasis = 4, norder = 4) {

  # Check if X is a 3-dimensional array
  if (length(dim(X)) != 3) {
    stop("X must be a 3-dimensional array, where the first dimension
         represent the number of observations, the second dimension
         represent the number of time points, and the third dimension
         reepresent the number of predictor variables.")
  }

  # Check if y is a vector of appropriate length
  if (!is.vector(y) | length(y) != dim(X)[1]) {
    stop("y must be a vector of length equal to n, the number of
    observations in X.")
  }

  # Check if H is a positive integer
  if (H <= 0 | H != floor(H)) {
    stop("H must be a positive integer.")
  }

  # Check if H is one number
  if (length(H) > 1) {
    stop("H must be one number.")
  }

  # Check if H is greater than 1 and less than the length of y
  if (H == 1 | H >= length(y)) {
    stop("H must be an integer that is at least 2 and less than the
         length of y.")
  }

  # Check if nbasis is a positive integer
  if (nbasis <= 0 | nbasis != floor(nbasis)) {
    stop("nbasis must be a positive integer.")
  }

  # Check if nbasis >= norder
  if (norder > nbasis) {
    stop(paste("nbasis must be at least norder; nbasis = ",nbasis, ",
               norder = ",norder, ""))
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
  databasis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis,
                                         norder = norder)

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

  # Compute the proportion of observations and means within each slice
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

#' Generate Functional data using pre-specified coefficients
#'
#' \code{fundata} generates functional data using either B-spline or Fourier
#' basis functions.
#'
#' This function generates functional data by expanding a set of coefficients
#' (`eta`) onto a chosen basis (either B-spline or Fourier).  Then, functional
#' principal component analysis (FPCA) is applied to obtain a low-dimensional
#' representation of the data.
#'
#' @param n An integer specifying the number of observations, i.e., sample size.
#' @param p An integer specifying the number of functional predictors.
#' @param nbasis An integer specifying the number of basis functions.
#' @param tt A numeric vector of time points at which the functional data is
#'     evaluated.
#' @param basisname A character string for the type of basis functions to use,
#'     either 'bspline' or 'fourier'.  Default is 'bspline'.
#' @param eta A matrix of coefficients of dimension \code{nbasis * n * p}.
#' @param norder An integer specifying the order of B-splines, which is one
#'     higher than their degree.  The default of 4 gives cubic splines.
#'
#' @return \code{fundata} returns a list containing:
#'    \itemize{
#'        \item \code{X}: An array of generated functional data of dimension
#'            \code{n * nt * p}, where \code{nt} denotes the number of time
#'            points.
#'        \item \code{Xc}: An array of centered functional data of dimension
#'            \code{n * nt * p}, where \code{nt} denotes the number of time
#'            points.
#'        \item \code{mfpca.scores}: A matrix of FPCA scores.
#'.       \item \code{xcoefs}: The array of the coefficients used for generating
#'        the functional data of dimension \code{nbasis * n * p}.  Note that, if
#'        Fourier basis is used and `nbasis` is even number, then it is rounded
#'        up to the nearest odd integer to preserve the pairing of sine and cosine
#'        functions.
#'        \item \code{xcoefs.mat}: The matrix version of the coefficients with
#'        dimension \code{n * (p * nbasis)}.  Note that, if Fourier basis is used
#'        and `nbasis` is even number, then it is rounded up to the nearest odd
#'        integer to preserve the pairing of sine and cosine functions.
#'        \item \code{basis}: The basis object used for functional data generation.
#'    }
#'
#' @examples
#' # Example 1: B-spline
#' # set the parameters
#' n <- 100
#' p <- 3
#' nbasis <- 4
#' tt <- seq(0, 1, length.out = 100)
#' eta.mat <- mvtnorm::rmvnorm(n, mean = rep(0, p * nbasis))
#' eta <- array(eta.mat, dim = c(nbasis, n, p))
#' # create the functional predictors
#' data <- fundata(n, p, nbasis, tt, 'bspline', eta)
#' str(data)
#' X <- data$X
#' # plot the first functional predictor for illustration
#'  fda::matplot(tt, t(X[, , 1]), type = "l", lty = 1, col = 1:n,
#'  xlab = "Time", ylab = "Value", main = paste("Functional Predictor", 1))
#'
#' # Example 2: Fourier
#' # set the parameters
#' n <- 100
#' p <- 5
#' nbasis <- 4
#' tt <- seq(0, 1, length.out = 100)
#' SigmaCov <- matrix(0, nrow = p * nbasis, ncol = p * nbasis)
#' for (j in 1:p) {
#'  index.j <-(((j - 1) * nbasis + 1):(j * nbasis))
#'  diag(SigmaCov[index.j, index.j]) <- c(2, 1, 1/2, 1/4)
#' }
#' eta.mat <- mvtnorm::rmvnorm(n, mean = rep(0, p * nbasis), sigma = SigmaCov)
#' eta <- array(eta.mat, dim = c(nbasis, n, p))
#' # create the functional predictors
#' data <- fundata(n, p, nbasis, tt, 'fourier', eta)
#' X <- data$X
#' # plot the first functional predictor for illustration
#' fda::matplot(tt, t(X[, , 1]), type = "l", lty = 1, col = 1:n,
#'     xlab = "Time", ylab = "Value", main = paste("Functional Predictor", 1))
#'
#' @export
fundata <- function(n, p, nbasis, tt, basisname = 'bspline', eta, norder = 4) {

  # Check if n is a single integer number
  if (length(n) != 1 | n != round(n) | n <= 0) {
    stop("Parameter 'n' must be a single positive integer number.")
  }

  # Check that p is a single number
  if (length(p) != 1) {
    stop("Parameter 'p' must be a single number.")
  }

  # Check that p is a positive integer
  if (p != round(p) | p <=0) {
    stop("Parameter 'p' must be a positive integer number.")
  }

  # Check that nbasis is a single number
  if (length(nbasis) != 1) {
    stop("Parameter 'nbasis' must be a single number.")
  }

  # Check that nbasis is a positive integer
  if (nbasis != round(nbasis) | nbasis <=0) {
    stop("Parameter 'nbasis' must be a positive integer number.")
  }

  # Check if norder if smaller than nbaasis - for Bsplines
  if (basisname == 'bspline' && norder > nbasis) {
    stop(paste("nbasis must be at least norder; nbasis = ",nbasis, ",
               norder = ",norder, ""))
  }

  # Check if time is a vector
  if (!is.vector(tt) | length(tt) == 1) {
    stop("time is a vector of length more than 1.")
  }

  # Check if the dimension of eta is n x p x nbasis
  if (!is.array(eta) | dim(eta)[1] != nbasis | dim(eta)[2] != n |
      dim(eta)[3] != p) {
    stop("Parameter 'eta' must be a numeric array with dimensions
         n x p x nbasis.")
  }

  # Step 1: Initialize the functional data array
  nt <- length(tt)
  X <- array(0, dim = c(n, nt, p))
  Xc <- array(0, dim = c(n, nt, p))

  # Step 2: Set up the basis functions
  if (basisname == 'bspline') {
    basis <- fda::create.bspline.basis(c(0, 1), norder = norder, nbasis = nbasis)
    st <- as.matrix(fda::eval.basis(tt, basis))
  } else if (basisname == 'fourier') {
    basis <- fda::create.fourier.basis(c(0, 1), nbasis = nbasis)
    nbasis_old <- nbasis
    nbasis <- basis$nbasis
    st <- as.matrix(fda::eval.basis(tt, basis))
  }

  # Step 3: Set up coefficients
  if (basisname == 'fourier') {
    if (nbasis_old %% 2 == 0) {
      xcoefs_extended <- array(0, c(nbasis, n, p))
      xcoefs_extended[2:nbasis, , ] <- eta
      eta <- xcoefs_extended
    }
  }

  # Step 4: Transform eta into a n * (p * nbasis) matrix
  eta.mat <- numeric()
  for (k in 1:p) {
    eta.mat <- cbind(eta.mat, t(eta[, , k]))
  }

  # Step 5: Generate functional data for each predictor
  for(j in 1:p) {
    index.j <- (((j - 1) * nbasis + 1):(j * nbasis))
    X[, , j] <- eta.mat[, index.j] %*% t(st)
    Xc[, , j] <- X[, , j] - matrix(rep(apply(X[, , j], 2, mean), n),
                                   nrow = n, byrow = T)
  }

  # Step 6: Perform Functional Principal Component Analysis (FPCA)
  pca.out <- fpca(list(coef = eta, basis = basis), basisname)
  mfpca.scores <- pca.out$pred

  # Return output as a list
  out <- list(X = X, Xc = Xc, mfpca.scores = mfpca.scores, xcoefs = eta,
              xcoefs.mat = eta.mat, basis = basis)
}

#' Generate Functional data
#'
#' \code{genfundata} generates functional data using either B-spline or Fourier
#' basis functions.
#'
#' This function constructs functional data by expanding a set of coefficients
#' onto a chosen basis (either B-spline or Fourier).  Then, functional principal
#' component analysis (FPCA) is applied to obtain a low-dimensional represetation
#' of the data.  If the coefficients are not provided, then they are randomly
#' generated.
#'
#' @param n Integer. The number of observations, i.e., sample size.
#' @param p Integer. The number of functional predictors.
#' @param nbasis Integer. The number of basis functions.
#' @param tt A numeric vector of time points at which the functional data is
#'     evaluated.
#' @param basisname A character string for the type of basis functions to use,
#'     either 'bspline' or 'fourier'.  Default is 'bspline'.
#' @param eta An optional array of coefficients of dimension \code{nbasis * n * p}.
#'     If not provided, the coefficients will be randomly generated.
#'
#' @return \code{genfundata} returns a list containing:
#'    \itemize{
#'        \item \code{X}: An array of generated functional data of dimension
#'            \code{n * nt * p}, where \code{nt} denotes the number of time
#'            points.
#'        \item \code{Xc}: An array of centered functional data of dimension
#'            \code{n * nt * p}, where \code{nt} denotes the number of time
#'            points.
#'        \item \code{mfpca.scores}: A matrix of FPCA scores.
#'.       \item \code{xcoefs}: An array of the coefficients used for generating
#'        the functional data of dimension \code{nbasis * n * p}.
#'        \item \code{xcoefs.mat}: A matrix version of the coefficients with
#'        dimension \code{n * (p * nbasis)}.
#'        \item \code{basis}: The basis object used for functional data generation.
#'    }
#'
#' @examples
#' # Example 1: Randomly generated coefficients
#' n <- 50
#' p <- 2
#' nbasis <- 5
#' tt <- seq(0, 1, length.out = 100)
#' data <- genfundata(n, p, nbasis, tt, 'bspline')
#' str(data)
#'
#' # Example 2: Coefficients are provided by the user
#' n <- 100
#' p <- 3
#' nbasis <- 4
#' tt <- seq(0, 1, length.out = 100)
#' eta.mat <- mvtnorm::rmvnorm(n, mean = rep(0, p * nbasis))
#' eta <- array(eta.mat, dim = c(nbasis, n, p))
#' data <- genfundata(n, p, nbasis, tt, 'bspline', eta = eta)
#' str(data)
#'
#' @export
genfundata <- function(n, p, nbasis, tt, basisname = 'bspline', eta = NULL) {

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

  # Check if time is a vector
  if (!is.vector(tt) | length(tt) == 1) {
    stop("time is a vector of length more than 1.")
  }

  # Check if the basis is Bspline of Fourier
  if (basisname != 'bspline' && basisname != 'fourier') {
    stop(paste('The basis type needs to be Bspline or Fourier.'))
  }

  if (!is.null(eta)) {
    if (!is.array(eta) | dim(eta)[1] != nbasis | dim(eta)[2] != n |
        dim(eta)[3] != p) {
      stop("Parameter 'eta' must be a numeric array with dimensions
         n x p x nbasis.")
    }
  }

  # Step 1: Initialize the functional data array
  nt <- length(tt)
  X <- array(0, c(n, nt, p))
  Xc <- array(0, c(n, nt, p))

  # Step 2: Set up the basis functions
  if (basisname == 'bspline') {
    basis <- fda::create.bspline.basis(c(0, 1), nbasis = nbasis)
  } else if (basisname == 'fourier') {
    basis <- fda::create.fourier.basis(c(0, 1), nbasis = nbasis)
    nbasis_old <- nbasis
    nbasis <- basis$nbasis
  }

  # If coefficients for the functional predictors are not provided
  if (is.null(eta)) {
  # Step 3: Generate random coefficients for functional predictors
  xcoefs <- array(0, c(n, nbasis, p))
  if (basisname == 'fourier') {
    xcoefs[, 1, ] <- rnorm(n * p, 0, 2)
    xcoefs[, 2, ] <- rnorm(n * p, 0, 1 / 4)
    xcoefs[, 3, ] <- rnorm(n * p, 0, 1 / 2)
    xcoefs[, 4, ] <- rnorm(n * p, 0, 1 / 4)
    if (nbasis >= 5) {
      for (i in 5:nbasis){
        xcoefs[, i, ] <- rnorm(n * p, 0, 1 / 2)
      }
    }
  } else if (basisname == 'bspline') {
    xcoefs[, 1, ] <- runif(n * p, -4, 4)
    xcoefs[, 2, ] <- runif(n * p, -1, 1)
    xcoefs[, 3, ] <- runif(n * p, -1 / 2, 1 / 2)
    xcoefs[, 4, ] <- runif(n * p, -1 / 2, 1 / 2)
    if (nbasis >= 5) {
      for (j in 5:nbasis) {
        xcoefs[, j, ] <- runif(n * p, -2, 2)
      }
    }

    # Normalize B-spline coefficients
    for (i in 1:p) {
      Xfd <- fda::fd(t(xcoefs[, , i]), basis)
      fd.sd <- sd(diag(fda::inprod(Xfd, Xfd)))
      xcoefs[, , i] <- xcoefs[, , i] / fd.sd
    }
  }

  # Step 4: Orthogonalize coefficients using Q-matrix
  Q <- qmat(n)
  zcoefs <- xcoefs
  for (i in 1:p) {
    xcoefs[, , i] <- Q %*% zcoefs[, , i]
  }
  xcoefs <- aperm(xcoefs, c(2, 1, 3))
  } else {
    xcoefs <- eta
    if (nbasis_old %% 2 == 0) {
      xcoefs_extended <- array(0, c(nbasis, n, p))
      xcoefs_extended[2:nbasis, , ] <- eta
      xcoefs <- xcoefs_extended
    }
  }

  # Step 5: Perform Functional Principal Component Analysis (FPCA)
  pca.out <- fpca(list(coef = xcoefs, basis = basis), basisname)
  mfpca.scores <- pca.out$pred

  # Step 6: Generate functional data for each predictor
  for (i in 1:p) {
    Xfd <- fda::fd(xcoefs[, , i], basis)
    X[, , i] <- t(fda::eval.fd(tt, Xfd))
    Xc[, , i] <- X[, , i] - matrix(rep(apply(X[, , i], 2, mean), n),
                                   nrow = n, byrow = T)
  }

  if (!is.null(eta) & nbasis_old %% 2 == 0) {
    xcoefs <- xcoefs[-1, , ]
  }

  # Create the xcoefs as a n x (p * nbasis) matrix
  xcoefs.mat <- numeric()
  for (k in 1:p) {
    xcoefs.mat <- cbind(xcoefs.mat, t(xcoefs[, , k]))
  }

  # Return output as a list
  out <- list(X = X, Xc = Xc, mfpca.scores = mfpca.scores, xcoefs = xcoefs,
              xcoefs.mat = xcoefs.mat, basis = basis)
  return(out)
}


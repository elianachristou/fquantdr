#' Functional data
#'
#' \code{genfundata} generates functional data using either B-spline or Fourier
#' basis functions.
#'
#' This function constructs functional data using basis expansions.  The
#' coefficients for the basis functions are randomly generated, and an
#' orthogonalization step is applied using a Q-matrix.  Functional principal
#' component analysis (FPCA) is then performed to obtain a low-dimensional
#' representation of the data.
#'
#' @param n The number of observations, i.e., sample size.
#' @param p The number of functional predictors.
#' @param nbasis The number of basis functions.
#' @param tt A vector of time points at which the functional data is
#'     evaluated.
#' @param basisname A character string for the type of basis functions to use,
#'     either 'bspline' or 'fourier'.  Default is 'bspline'.
#'
#' @return \code{genfundata} returns a list containing:
#'    \itemize{
#'        \item \code{X}: An array of generated functional data of dimensions
#'            \code{n * nt * p} array, where \code{nt} denotes the number
#'            of time points.
#'        \item \code{mfpca.scores}: Matrix of FPCA scores.
#'    }
#'
#' @examples
#' n <- 50
#' p <- 2
#' nbasis <- 5
#' tt <- seq(0, 1, length.out = 100)
#' data <- genfundata(n, p, nbasis, tt, 'bspline')
#' str(data)
#'
#' @export
genfundata <- function(n, p, nbasis, tt, basisname = 'bspline') {
  # Step 1: Initialize the functional data array
  nt <- length(tt)
  X <- array(0, c(n, nt, p))

  # Step 2: Set up the basis functions
  if (basisname == 'bspline') {
    basis <- fda::create.bspline.basis(c(0, 1), nbasis = nbasis)
  } else if (basisname == 'fourier') {
    basis <- fda::create.fourier.basis(c(0, 1), nbasis = nbasis)
  }

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

  # Step 5: Perform Functional Principal Component Analysis (FPCA)
  pca.out <- fpca(list(coef = xcoefs, basis = basis), basisname)
  mfpca.scores <- pca.out$pred

  # Step 6: Generate functional data for each predictor
  for (i in 1:p) {
    Xfd <- fda::fd(xcoefs[, , i], basis)
    X[, , i] <- t(fda::eval.fd(tt, Xfd))
  }

  # Return output as a list
  out <- list(X = X, mfpca.scores = mfpca.scores, xcoefs = xcoefs)
  return(out)
}


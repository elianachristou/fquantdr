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
#' @param Xc A \code{n x t x p} array, where n is the number of observations,
#'  t is the number of time points, and p is the number of variables.
#' @param y A numeric vector of length \code{n} representing the response
#'     variable.
#' @param H An integer representing the number of slices for the response
#'     variable.
#' @param nbasis An integer representing the number of basis functions for
#'     the B-spline basis.
#'
#' @return \code{mfsir} computes the new sufficient predictors and returns:
#' \item{phi}{The eigenvectors of the transformation matrix.}
#' \item{betas}{The regression coefficients in the B-spline basis.}
#' \item{eigvalues}{The eigenvalues of the transformation matrix.}
#' \item{xfd.coef}{The coordinates of the smoothed data.}
#' \item{gx}{The block diagonal Gram matrix.}
#'
#' @examples
#' # set the parameters
#' n <- 100
#' p <- 5
#' nbasis <- 4
#' H <- 10
#' time_points <- seq(0, 1, length.out = 101)
#' eta <- matrix(stats::rnorm(n * p * q), nrow = n, ncol = p * q)
#'
#' # Generate functional data
#' result <- fundata(n, p, q, time_points, eta)
#'
#' # Centered functional predictors
#' Xc <- result$cg
#'
#' P <- eigen(cov(eta))$vectors
#' Q <- diag(eigen(cov(eta))$values)
#' # This is the inner products <b1, X>, <b2, X>, ...
#' mfpca.scores <- eta %*% P
#' # Define true inner products for functional central subspace
#' uu_fcs <- mfpca.scores[, 1:K]
#' # Generate error for data
#' error <- rnorm(n)
#' # Generate response from functional data and error
#' y <- 3 * mfpca.scores[, 1] + error
#'
#' # Run mfsir
#' result <- mfsir(Xc, y, H, nbasis)
#'
#' @noRd
mfsir <- function(Xc, y, H, nbasis) {
  # Check if Xc is a 3-dimensional array
  if (length(dim(Xc)) != 3) {
    stop("Xc must be a 3-dimensional array")
  }

  # Check if y is a vector of appropriate length
  if (!is.vector(y) || length(y) != dim(Xc)[1]) {
    stop("y must be a vector of length equal to the number of
         observations in Xc")
  }

  # Check if H is an integer
  if (!is.numeric(H) || H <= 0 || H != floor(H)) {
    stop("H must be a positive integer")
  }

  # Check if nbasis is an integer
  if (!is.numeric(nbasis) || nbasis <= 0 || nbasis != floor(nbasis)) {
    stop("nbasis must be a positive integer")
  }

  n <- dim(Xc)[1] # Number of observations
  p <- dim(Xc)[3] # Number of variables

  # Creates a B-spline basis for smoothing
  databasis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)

  # Check if databasis is of class basisfd
  if (!inherits(databasis, "basisfd")) {
    stop("databasis must be of class basisfd")
  }

  # Initialize an empty vector to store coefficients
  xfd.coef <- numeric() # retrieves the coordinates as an n x (p * q) matrix

  # Loop over each variable to smooth and center the data
  for (k in 1:p) {
    xfdk <- fda::smooth.basis(seq(0, 1, length.out = dim(Xc)[2]), t(Xc[, , k]),
                              databasis)$fd
    xfdk <- fda::center.fd(xfdk)
    xk.coef <- t(xfdk$coef)
    xfd.coef <- cbind(xfd.coef, xk.coef)
  }

  #compute block diagonal gram matrix gx
  gi <- gramatrix(nbasis, databasis)
  gx <- matrix(0, nrow = p * (nbasis), ncol = p * (nbasis))

  # Fill the block diagonal Gram matrix
  for (i in 1:p) {
    index <- ((i - 1) * (nbasis) + 1):(i * (nbasis))
    gx[index, index] <- gi
  }

  # Compute the covariance matrix of the smoothed data coefficients
  sigmaxx <- (t(xfd.coef) %*% xfd.coef) / n

  # Compute the square root of the Gram matrix
  gx.sqrt <- mppower(gx, 1/2, 10e-7)

  # Compute the product of gx.sqrt, sigmaxx, and gx.sqrt
  store <- gx.sqrt %*% sigmaxx %*% gx.sqrt

  # Compute the inverse square root of the product matrix
  store.inv.sqrt <- rigpower(store, -1/2, 0.05)

  # Compute the mean of the coefficients
  mucoef <- apply(xfd.coef, 2, mean)

  # Discretize the response variable into H slices
  ydis <- discretize (y, H)

  # Compute the slice probabilities
  prob <- slprob(ydis, H)

  # Compute the average slice values
  avg.slav <- slav(xfd.coef, ydis, H)

  # Initialize Lambda1 matrix
  Lambda1 <- matrix(0, p * nbasis, p * nbasis)

  # Compute Lambda1 using the slice probabilities and average slice values
  for(i in 1:H) {
    Lambda1 <- Lambda1 + prob[i] * avg.slav[i, ] %*% t(avg.slav[i, ])}

  # Compute the cross-covariance matrix
  sigmaxxy <- Lambda1 - mucoef %*% (t(mucoef))

  # Compute the transformation matrix
  mat <- store.inv.sqrt %*% gx.sqrt

  # Compute the matrix to be decomposed
  out <- mat %*% sigmaxxy %*% t(mat)

  # Ensure the matrix is symmetric
  out <- symmetry(out)

  # Perform eigen decomposition
  eigendecom <- eigen(out)

  # Extract eigenvalues and eigenvectors
  eigvalues <- eigendecom$values
  phi <- eigendecom$vec

  # Compute the inverse square root of the Gram matrix
  gx.inv.sqrt <- rigpower(gx, -1/2, 0.05)

  # Compute the regression coefficients
  betas <- gx.inv.sqrt %*% store.inv.sqrt %*% phi

  # Return the results as a list
  return(list(phi = phi, betas = betas, eigvalues = eigvalues, xfd.coef = xfd.coef, gx = gx))
}

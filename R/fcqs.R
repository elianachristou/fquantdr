#' Functional Central Quantile Subspace
#'
#' \code{fcqs} estimates the directions of the functional central quantile
#' subspace, extending the central quantile subspace to functional data.
#'
#' This function computes the directions that span the \eqn{\tau}th functional
#' central quantile subspace. These directions represent functions that can
#' be linearly applied via the inner product to given predictors to reduce the
#' dimension of infinitely-dimensional functional data without losing any
#' information required for accurate quantile regression of the functional
#' data.
#'
#' @param x A 3-dimensional array (\code{n x nt x p}), where n is the number
#'     of observations, nt is the number of time points, and p is the number
#'     of predictor variables.
#' @param y A numeric vector of length \code{n} representing the response
#'     variable.
#' @param time A numeric vector of length \code{nt} of time points at which
#'     the functional data is evaluated.
#' @param nbasis The number of basis functions for smoothing functional data.
#' @param tau A quantile level, a number strictly between 0 and 1.
#' @param d_tau The dimension of the functional central quantile subspace.
#'      It should be an integer between 1 and `p`.
#'
#' @return `fcqs` computes the directions of the functional central quantile
#'      subspace (FCQS) and returns:
#'      \item{betacoef}{The functional parameters that span the FCQS}
#'      \item{yhat}{The estimated responses, resulting as the inner product
#'      between `betacoef` and `x`.}
#'
#' @export
#'
#' @examples
#' # Parameters for generating functional data
#' n <- 100
#' p <- 5
#' q <- 4
#' time <- seq(0, 1, length.out = 101)
#' eta <- matrix(stats::rnorm(n * p * q), nrow = n, ncol = p * q)
#'
#' # Generate functional data
#' result <- fundata(n, p, q, time, eta)
#'
#' # Centered functional predictors
#' xc <- result$xc
#'
#' # Further parameters for FCQS
#' H <- 10
#' tau <- 0.1
#' K <- 1
#' d_tau <- 1
#' nbasis <- q
#' d_DR <- 1
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
#' # Run FCQS function
#' fcqs(xc, y, time_points, q, nbasis, tau, d_tau, H, d_DR)

fcqs <- function(x, y, time, nbasis, tau = 0.5, d_tau) {

  # Check if y is a univariate response
  if (!is.vector(y)) {
    stop("y should be a vector representing a univariate response.")
  }

  # Check if x is an array
  if (!is.array(x)) {
    stop(paste('x needs to be an array; convert it to a (n x nt x p) array,
               where n is the sample size, nt is the number of time points,
               and p is the number of predictors.'))
  }

  # Check if x is a 3-dimensional array
  if (length(dim(X)) != 3) {
    stop("X must be a 3-dimensional array, where the first dimension
         represent the number of observations, the second dimension
         represent the number of time points, and the third dimension
         reepresent the number of predictor variables.")
  }

  # Check if time is a univariate vector
  if (!is.vector(time)) {
    stop("The input 'time' should be a vector.")
  }

  # Check that the dimensions agree
  if(dim(x)[1] != length(y)) {
    stop("y and x should have the same number of observations.")
  }

  if (dim(x)[2] != length(time)) {
    stop(paste("x and time should have the same number of time points."))
  }

  # Check if tau is between 0 and 1
  if(tau < 0 | tau > 1) {
    stop("The quantile level, tau, must be strictly between 0 and 1.")
  }

  # Set the parameters
  n <- dim(x)[1]
  nt <- dim(x)[2]
  p <- dim(x)[3]
  H <- max(10, 2 * p / n)

  # Create basis for functional data using splines
  databasis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
  # Initialize array of predictor coordinates (n x (p * q) dimension)
  xcoef <- numeric()
  # Retrieves the coefficients as a (q x n x p) array for the fRegress function
  xcoef_array <- array(0, c(nbasis, n, p))
  # Loop for all predictors
  for (k in 1:p) {
    # Smooth functional predictor using given basis
    xfdk <- fda::smooth.basis(time_points, t(x[, , k]), databasis)$fd
    # Center data
    xfdk <- fda::center.fd(xfdk)
    # Get coefficients
    xk.coef <- t(xfdk$coefs)
    # Create array from coefficients
    xcoef_array[, , k] <- t(xk.coef)
    # Add coefficients to xcoef matrix
    xcoef <- cbind(xcoef, xk.coef)
  }

  # Compute block diagonal gram matrix gx
  gi <- gramatrix(nbasis, databasis)
  # Initialize gx as (p * nbasis) x (p * nbasis) square matrix
  gx <- matrix(0, nrow = p * (nbasis), ncol = p * (nbasis))
  # Loop through predictors
  for (i in 1:p) {
    # Get indices for (nbasis x nbasis) block matrix for predictor in gx
    index <- ((i - 1) * (nbasis) + 1):(i * (nbasis))
    # Add block matrix at diagonal index in gx
    gx[index, index] <- gi
  }

  # Run Functional Sliced Inverse Regression
  out.mfsir <- mfsir(x, y, H, nbasis)
  beta.mfsir <- out.mfsir$betas[, 1:d_DR]
  vv <- xcoef %*% gx %*% beta.mfsir

  # Run FCQS methodology to first find the initial vector
  # Estimate the conditional quantile
  red_dim <- floor(0.1 * n) # Find length of 10% of observations
  # Get middle 80% of response
  index_y <- order(y)[red_dim:(n - red_dim)]
  # Get initial bandwidth
  h <- KernSmooth::dpill(vv[index_y, ], y[index_y])
  # Transform bandwidth for given quantile
  h <- 4 * h * (tau * (1 - tau) / (stats::dnorm(stats::qnorm(tau))) ^ 2) ^ 0.2
  # Use alternative method for bandwidth calculation on error
  if (h == 'NaN') {
    h <- 1.25 * max(n ^ (-1 / (d_tau + 4)), min(2, stats::sd(y)) *
                    n ^ (-1 / (d_tau + 4)))
  }
  # Scale bandwidth by 3
  h <- 3 * h
  # Get qhat through LLQR function
  qhat <- quantdr::llqr(vv, y, tau = tau, h = h)$ll_est

  # Fit a simple linear regression from qhat on xfd
  # Create functional data object
  x.fd <- fda::fd(xcoef_array, databasis)
  fcqs.out <- sonf(qhat, x.fd, dev2_penalty = TRUE, lambda = 1e-4)
  # Compute initial inner product
  initial_inner_prod <- fcqs.out$yhat
  # Compute initial beta coefficients
  initial_beta_coef <- fcqs.out$beta_coef

  # Initialize beta coefficients
  beta_vector <- initial_beta_coef
  # Initialize inner product values
  inner_prod <- initial_inner_prod
  # Construct more vectors
  for (j in 1:(min(p * q, 40) - 1)) {
    # Construct new beta vector using the inner product and x coefficients
    new_beta <- apply(matrix(rep(as.vector(inner_prod), 2), n, p * q) * xcoef,
                      2, mean)
    # Combine new beta vector and initial beta vector
    beta_vector <- cbind(beta_vector, new_beta)

    # Need to create the new inner_prod
    inner_prod <- xcoef %*% gx %*% new_beta
  }

  # Eigenfunction decomposition
  # Compute Moore-Penrose inverse of gx to the .5 power
  gx.half = mppower(gx, 0.5, 1e-8)
  # Compute Moore-Penrose inverse of gx to the -.5 power
  gx.inv.half = mppower(gx, -0.5, 1e-8)
  # Compute matrix A to find eigenvectors of M
  A =  mppower(1 / n * gx.half %*% t(xcoef) %*% xcoef %*% gx.half, -0.5, 1e-8)
  M <- A %*% gx.half %*% beta_vector %*% t(beta_vector) %*% gx.half %*% A
  # Obtain eigenvectors and eigenvalues of symmetric M
  gg <- eigen(M, symmetric = T)$vectors
  # Compute directions of central quantile subspace
  vv <- gx.inv.half %*% A %*% gg[, 1:d_tau]
  vv2 <- xcoef %*% gx.half %*% A %*% gg[, 1:d_tau]

  # Return directions of FCQS
  list(ffun = vv, yhat = vv2)
}

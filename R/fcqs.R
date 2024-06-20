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
#' @param Xc An n x t x p array of functional data where n represents the
#' number of observations, t represents the number of time points, and p is the
#' number of functional predictors
#' @param y A vector representing the univariate response
#' @param time_points A vector of the time points spanned by the functional
#' predictors
#' @param q A number used to define the dimension during KL expansion
#' @param nbasis The number of basis functions for smoothing functional data
#' @param tau The quantile level, between 0 and 1, which we are calculating
#' the FCQS for
#' @param d_tau The dimension of the functional central quantile subspace (FCQS)
#' @param H The number of slices for functional sliced inverse regression
#' @param d_DR The dimension of the functional central subspace (FCS)
#'
#' @return The directions of the functional central quantile subspace, which can
#' be used to sufficiently reduce the dimension of the functional predictors for
#' quantile regression
#' @export
#'
#' @examples
fcqs <- function(Xc, y, time_points, q, nbasis, tau, d_tau, H, d_DR) {

  if (!is.array(Xc)) {
    stop(paste('Xc needs to be an array; convert it to a n x nt x p array,',
               'where n is the sample size, nt is the number of time points,',
               'and p is the number of predictors'))
  }

  if (!is.vector(y)) {
    stop("y should be a vector representing a univariate response.")
  }

  if (!is.vector(time_points)) {
    stop("time_points should be a vector.")
  }

  if (dim(Xc)[2] != length(time_points)) {
    stop(paste('Xc needs to be an n x nt x p array, where n is the sample',
               'size, nt is the number of time points, and p is the number of',
               'predictors'))
  }

  if(dim(Xc)[1] != length(y)) {
    stop("y and Xc should show the same number of observations.")
  }

  if(tau < 0 | tau > 1) {
    stop("The quantile level, tau, must be between 0 and 1.")
  }

  # Get dimensions of predictors
  # n = number of observations
  # ntx = number of time points
  # p = number of predictors
  n <- dim(Xc)[1]
  ntx <- dim(Xc)[2]
  p <- dim(Xc)[3]

  # Create basis for functional data using splines
  databasis <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
  # Initialize array of predictor coordinates (n x (p * q) dimension)
  xcoef <- numeric()
  # Retrieves the coefficients as a (q x n x p) array for the fRegress function
  xcoef_array <- array(0, c(nbasis, n, p))
  # Loop for all predictors
  for (k in 1:p) {
    # Smooth functional predictor using given basis
    xfdk <- fda::smooth.basis(time_points, t(Xc[, , k]), databasis)$fd
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
  out.mfsir <- mfsir(Xc, y, H, nbasis)
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
  h <- 4 * h * (tau * (1 - tau) / (dnorm(qnorm(tau))) ^ 2) ^ 0.2
  # Use alternative method for bandwidth calculation on error
  if (h == 'NaN') {
    h <- 1.25 * max(n ^ (-1 / (d_tau + 4)), min(2, sd(y)) *
                    n ^ (-1 / (d_tau + 4)))
  }
  # Scale bandwidth by 3
  h <- 3 * h
  # Get qhat through LLQR function
  qhat <- quantdr::llqr(vv, y, tau = tau, h = h)$ll_est

  # Fit a simple linear regression from qhat on xfd
  # Create functional data object
  x.fd <- fd(xcoef_array, databasis)
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
  # Compute Moor-Penrose inverse of gx to the -.5 power
  gx.inv.half = mppower(gx, -0.5, 1e-8)
  # Compute matrix A to find eigenvectors of M
  A =  mppower(1 / n * gx.half %*% t(xcoef) %*% xcoef %*% gx.half, -0.5, 1e-8)
  M <- A %*% gx.half %*% beta_vector %*% t(beta_vector) %*% gx.half %*% A
  # Obtain eigenvectors and eigenvalues of symmetric M
  gg <- eigen(M, sym = T)$vectors
  # Compute directions of central quantile subspace
  vv <- gx.inv.half %*% A %*% gg[, 1:d_tau]
  vv2 <- xcoef %*% gx.half %*% A %*% gg[, 1:d_tau]

  # Return directions of FCQS
  list(vv = vv, vv2 = vv2)
}

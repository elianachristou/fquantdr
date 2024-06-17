fcqs <- function(Xc, y, time_points, q, nbasis, tau, d_tau, H, d_DR){

  if (!is.array(Xc)) {
    stop(paste('Xc needs to be an array; convert it to a n x nt x p array,',
               'where n is the sample size, nt is the number of time points,',
               'and p is the number of predictors'))
  }

  if (dim(Xc)[2] != length(time_points)){
    stop(paste('Xc needs to be an n x nt x p array, where n is the sample size,',
               'nt is the number of time points, and p is the number of',
               'predictors'))
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
    # Get coefficients (Note: Should this be coefs instead of coef?)
    xk.coef <- t(xfdk$coef)
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

  # run fcqs to first find the initial vector
  # estimate the conditional quantile
  # find how many observations correspond to a 10%
  red_dim <- floor(0.1 * n)
  index_y <- order(y)[red_dim:(n - red_dim)]
  # estimates bandwidth h
  h <- KernSmooth::dpill(vv[index_y, ], y[index_y])
  # changes bandwidth calculation based on tau
  h <- 4 * h * (tau * (1 - tau) / (dnorm(qnorm(tau)))^2)^.2
  # deal with a non-numeric bandwidth values
  if (h == 'NaN') {
    h <- 1.25 * max(n^(-1 / (d_tau + 4)), min(2, sd(y)) * n^(- 1 / (d_tau + 4)))
  }
  h <- 3 * h
  # computes conditional quantile estimate
  qhat <- quantdr::llqr(vv, y, tau = tau, h = h)$ll_est

  # fit a simple linear regression from qhat on xfd
  # create functional data object
  x.fd <- fd(xcoef_array, databasis)
  fcqs.out <- sonf(qhat, x.fd, dev2_penalty = TRUE, lambda = 1e-4)
  # compute initial inner product
  initial_inner_prod <- fcqs.out$yhat
  # compute initial beta coefficients
  initial_beta_coef <- fcqs.out$beta_coef

  # initialize beta coefficients
  beta_vector <- initial_beta_coef
  # initialize innner product values
  inner_prod <- initial_inner_prod
  # construct more vectors
  for (j in 1:(min(p * q, 40) - 1)) {
    # construct new beta vector using the inner product and x coefficients
    new_beta <- apply(matrix(rep(as.vector(inner_prod), 2), n, p * q) * xcoef,
                      2, mean)
    # combine new beta vector and initial beta vector
    beta_vector <- cbind(beta_vector, new_beta)

    # need to create the new inner_prod
    inner_prod <- xcoef %*% gx %*% new_beta
  }

  # eigenfunction decomposition
  # compute Moore-Penrose inverse of gx to the .5 power
  gx.half = mppower(gx, 0.5, 10^(-8))
  # compute Moor-Penrose inverse of gx to the -.5 power
  gx.inv.half = mppower(gx, -0.5, 10^(-8))
  # Compute matrix A to find eigenvectors of M
  A =  mppower(1/n*gx.half %*% t(xcoef) %*% xcoef %*% gx.half,-0.5,10^(-8))
  M <- A %*% gx.half %*% beta_vector %*% t(beta_vector) %*% gx.half %*% A
  # obtain eigenvectors and eigenvalues of symmetric M
  gg <- eigen(M, sym = T)$vectors
  #
  vv <- gx.inv.half %*% A %*% gg[, 1:d_tau]
  vv2 <- xcoef %*% gx.half %*% A %*% gg[,1:d_tau]

  return(list(vv = vv, vv2 = vv2))
}

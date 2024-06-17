fcqs <- function(Xc, y, t, q, nbasis, tau, d_tau, H, d_DR){

  if (!is.array(Xc)) {
    stop(paste('Xc needs to be an array; convert it to a n x nt x p array,',
               'where n is the sample size, nt is the number of time points,',
               ' and p is the number of predictors'))
  }

  if (dim(Xc)[2] != length(t)){
    stop(paste('Xc needs to be an n x nt x p array, where n is the sample size,',
               ' nt is the number of times points, and p is the number of
               ', 'predictors'))
  }

  n <- dim(Xc)[1]
  ntx <- dim(Xc)[2]
  p <- dim(Xc)[3]

  databasis <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
  xcoef <- numeric() #retrieves the coodinates as an n x (p * q) matrix
  # retrieves the coefficients as a (q x n x p) array for the fRegress function
  xcoef_array <- array(0, c(nbasis, n, p))
  for (k in 1:p) {
    xfdk <- fda::smooth.basis(t, t(Xc[, , k]), databasis)$fd
    xfdk <- fda::center.fd(xfdk)
    xk.coef <- t(xfdk$coef)
    xcoef_array[, , k] <- t(xk.coef)
    xcoef <- cbind(xcoef, xk.coef)
  }

  #compute block diagonal gram matrix gx
  gi <- gramatrix(nbasis, databasis)

  gx <- matrix(0, nrow = p * (nbasis), ncol = p * (nbasis))
  for (i in 1:p) {
    index <- ((i - 1) * (nbasis) + 1):(i * (nbasis))
    gx[index, index] <- gi
  }

  # run fsir
  out.mfsir <- mfsir(Xc, y, H, nbasis)
  beta.mfsir <- out.mfsir$betas[, 1:d_DR]
  vv <- xcoef %*% gx %*% beta.mfsir

  # run fcqs to first find the initial vector
  # estimate the conditional quantile
  red_dim <- floor(0.1 * n) # find how many observations correspond to a 10%
  index_y <- order(y)[red_dim:(n - red_dim)]
  h <- KernSmooth::dpill(vv[index_y, ], y[index_y])
  h <- 4 * h * (tau * (1 - tau) / (dnorm(qnorm(tau)))^2)^.2
  if (h == 'NaN') {
    h <- 1.25 * max(n^(-1 / (d_tau + 4)), min(2, sd(y)) * n^(- 1 / (d_tau + 4)))
  }
  h <- 3 * h
  qhat <- llqr(vv, y, tau = tau, h = h)$ll_est

  # fit a simple linear regression from qhat on xfd
  x.fd <- fd(xcoef_array, databasis)
  fcqs.out <- sonf(qhat, x.fd, dev2_penalty = TRUE, lambda = 1e-4)
  initial_inner_prod <- fcqs.out$yhat
  initial_beta_coef <- fcqs.out$beta_coef # or betafd

  beta_vector <- initial_beta_coef
  inner_prod <- initial_inner_prod
  # construct more vectors
  for (j in 1:(min(p * q, 40) - 1)) {
    new_beta <- apply(matrix(rep(as.vector(inner_prod), 2), n, p * q) * xcoef,
                      2, mean)
    beta_vector <- cbind(beta_vector, new_beta)

    # need to create the new inner_prod
    inner_prod <- xcoef %*% gx %*% new_beta
  }

  # eigenfunction decomposition
  gx.half = mppower(gx, 0.5, 10^(-8))
  gx.inv.half = mppower(gx, -0.5, 10^(-8))
  A =  mppower(1/n*gx.half %*% t(xcoef) %*% xcoef %*% gx.half,-0.5,10^(-8))
  M <- A %*% gx.half %*% beta_vector %*% t(beta_vector) %*% gx.half %*% A
  gg <- eigen(M, sym = T)$vectors
  vv <- gx.inv.half %*% A %*% gg[, 1:d_tau]
  vv2 <- xcoef %*% gx.half %*% A %*% gg[,1:d_tau]

  return(list(vv = vv, vv2 = vv2))
}

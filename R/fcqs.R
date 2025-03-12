#' Functional Central Quantile Subspace
#'
#' \code{fcqs} estimates the directions of the functional central quantile
#' subspace.
#'
#' This function computes the directions that span the \eqn{\tau}th functional
#' central quantile subspace. These directions represent functions that can
#' be linearly applied via the inner product to replace the infinite-dimensional
#' functional predictors with a few finite predictors without losing important
#' information on the conditional quantiles while maintaining a flexible
#' nonparametric model. This methodology is introduced in Christou et al. (2025+).
#'
#' @param x A 3-dimensional array (\code{n x nt x p}), where n is the number
#'     of observations, nt is the number of time points, and p is the number
#'     of predictor variables.
#' @param y A numeric vector of length \code{n} representing the scalar
#'     response variable.
#' @param tt A numeric vector of length \code{nt} of time points at which
#'     the functional data is evaluated.
#' @param nbasis The number of basis functions for smoothing the functional
#'     data.
#' @param tau A quantile level, a number strictly between 0 and 1.
#' @param dtau The number of directions to extract, i.e., the dimension of
#'      the functional central quantile subspace.  It should be an integer
#'      between 1 and `p`.  If not provided, the function will return `p`
#'      directions
#'
#' @return `fcqs` computes the directions of the functional central quantile
#'      subspace (FCQS) and returns:
#'      \item{betacoef}{A matrix of size \code{nbasis x dtau}, representing
#'      the functional parameters that span the FCQS}
#'      \item{betax}{A numeric vector of length \code{n}, representing the
#'      sufficient predictor, calculated as the inner product between
#'      `betacoef` and `x`.}
#'
#' @references Christou, E., Solea, E., Wang, S., and Song, J. (2025+) Sufficient
#' Dimension Reduction for the Conditional Quantiles of Functional Data.
#' \emph{under revision}.
#'
#' @export
#'
#' @examples
#' # Set the parameters
#' n <- 100
#' p <- 5
#' nbasis <- 4
#' nt <- 100
#' tau <- 0.1
#' dtau <- 1
#' tt <- seq(0, 1, length.out = nt)
#' # Set the covariance matrix
#' SigmaCov <- matrix(0, p * nbasis, p * nbasis)
#' for (j in 1:p) {
#'  index.j <-(((j - 1) * nbasis + 1):(j * nbasis))
#'  diag(SigmaCov[index.j, index.j]) <- c(2, 1, 1/2, 1/4)
#' }
#' eta.mat <- mvtnorm::rmvnorm(n, mean = rep(0, p * nbasis), sigma = SigmaCov)
#' eta <- array(eta.mat, dim = c(nbasis, n, p))
#' # Generate functional data
#' data <- fundata(n, p, nbasis, tt, 'fourier', eta)
#' Xc <- data$Xc
#' # Generate y
#' error <- rnorm(n)
#' y <- 3 * data$mfpca.scores[, 1] + error
#' # Run FCQS function
#' output <- fcqs(Xc, y, tt, nbasis, tau, dtau)
#' # Compare the true and estimated predictors
#' mcorr(output$betax, mfpca.scores[, 1])

fcqs <- function(x, y, tt, nbasis, tau = 0.5, dtau = NULL) {

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
  if (length(dim(x)) != 3) {
    stop("x must be a 3-dimensional array, where the first dimension
         represent the number of observations, the second dimension
         represent the number of time points, and the third dimension
         reepresent the number of predictor variables.")
  }

  # Check if time is a univariate vector
  if (!is.vector(tt)) {
    stop("The input 'time' should be a vector.")
  }

  # Check that the dimensions agree
  if(dim(x)[1] != length(y)) {
    stop("y and x should have the same number of observations.")
  }
  if (dim(x)[2] != length(tt)) {
    stop(paste("x and time should have the same number of time points."))
  }

  # Check if tau is between 0 and 1
  if(tau <= 0 | tau >= 1) {
    stop("The quantile level, tau, must be strictly between 0 and 1.")
  }

  # Check that nbasis is a single number
  if (length(nbasis) != 1) {
    stop("The input 'nbasis' must be a single number.")
  }

  # Check that nbasis is a positive integer
  if (nbasis != round(nbasis) | nbasis <=0) {
    stop("The input 'nbasis' must be a positive integer number.")
  }

  # If dtau is missing, set it equal to p.  Otherwise, perform compatibility checks
  if (is.null(dtau)) {
    dtau <- dim(x)[3]
  } else {
   # Check if dtau is an integer
   if (floor(dtau) != dtau) {
     stop("The input 'dtau' must be an integer.")
   }

   # Check if dtau is between 1 and p
   if (dtau < 1 | dtau > dim(x)[3]) {
     stop("The input 'dtau' should be an integer between 1 and p, the
          number of predictor variables.")
   }
  }

  # Set the parameters
  n <- dim(x)[1]
  nt <- dim(x)[2]
  p <- dim(x)[3]
  H <- max(10, floor(2 * p / n))

  # Create basis for functional data using splines
  if (nbasis < 4) {
  databasis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis,
                                         norder = nbasis)
  } else {
    databasis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
  }

  # Retrieve the coefficients as a (nbasis x n x p) array for the sonf function
  xcoef <- numeric()
  xcoef_array <- array(0, c(nbasis, n, p))
  for (k in 1:p) {
    # Smooth functional predictor using given basis
    xfdk <- fda::smooth.basis(tt, t(x[, , k]), databasis)$fd
    # Center data and get coefficients
    xfdk <- fda::center.fd(xfdk)
    xk.coef <- t(xfdk$coefs)
    xcoef_array[, , k] <- t(xk.coef)
    xcoef <- cbind(xcoef, xk.coef)
  }

  # Compute block diagonal gram matrix gx
  gi <- gramatrix(nbasis, databasis)
  gx <- matrix(0, nrow = p * (nbasis), ncol = p * (nbasis))
  for (i in 1:p) {
    index <- ((i - 1) * (nbasis) + 1):(i * (nbasis))
    gx[index, index] <- gi
  }

  # Run Functional Sliced Inverse Regression
  out.mfsir <- mfsir(x, y, H, nbasis)
  beta.mfsir <- out.mfsir$betas[, 1:dtau]
  vv <- xcoef %*% gx %*% beta.mfsir

  # Run FCQS methodology to first find the initial vector
  ## Calculate the bandwidth and estimate the conditional quantile
  ## Trim the y for the bandwidth calculation
  red_dim <- floor(0.1 * n)
  index_y <- order(y)[red_dim:(n - red_dim)]
  h <- NULL
  h <- tryCatch({
    # Try the first method, using subset of pred and y
    h_value <- KernSmooth::dpill(vv[index_y, ], y[index_y])

    # Check if bandwidth is valid (not missig or non-positive)
    if (is.na(h_value) || (!missing(h_value) && h_value <= 0) || is.nan(h_value))
      stop("'bandwidth' must be strictly positive")

    h_value  # Return this value if no error
  }, error = function(e) {

    # Fallback to the alternative bandwidth calculation method using another tryCatch
    h_fallback <- tryCatch({
      # Attempt to calculate using the full dataset
      h_fallback <- KernSmooth::dpill(vv, y)

      # Check if fallback bandwidth is valid (not missing or non-positive)
      if (is.na(h_fallback) || (!missing(h_fallback) && h_fallback <= 0) || is.nan(h_fallback))
        stop("'fallback bandwidth' must be strictly positive")

      h_fallback  # Return fallback value if no error
    }, error = function(e2) {
      # If the fallback method also fails, use the final backup formula
      #message("Both primary and fallback bandwidth calculations failed. Using the backup formula.")

      # Final backup formula for bandwidth
      h <- 1.25 * max(n ^ (-1 / (dtau + 4)), min(2, stats::sd(y)) *
                        n ^ (-1 / (dtau + 4)))
      })

    h_fallback  # Return the fallback value or the backup formula value
  })

  h <- 3 * 4 * h * (tau * (1 - tau) / (stats::dnorm(stats::qnorm(tau))) ^ 2) ^ 0.2

  # Get qhat through LLQR function
  qhat <- quantdr::llqr(vv, y, tau = tau, h = h)$ll_est

  # Fit a simple linear regression from qhat on xfd
  x.fd <- fda::fd(xcoef_array, databasis)
  fcqs.out <- sonf(qhat, x.fd, dev2_penalty = TRUE, lambda = 1e-4)
  # Compute initial inner product
  initial_inner_prod <- fcqs.out$yhat
  # Compute initial beta coefficients
  initial_beta_coef <- fcqs.out$betacoef

  # Construct more vectors
  beta_vector <- initial_beta_coef
  inner_prod <- initial_inner_prod
  for (j in 1:(min(p * nbasis, 40) - 1)) {
    # Construct new beta vector using the inner product and x coefficients
    new_beta <- apply(matrix(rep(as.vector(inner_prod), p * nbasis), n,
                             p * nbasis) * xcoef, 2, mean)
    beta_vector <- cbind(beta_vector, new_beta)
    inner_prod <- xcoef %*% gx %*% new_beta
  }

  # Eigenfunction decomposition
  gx.half = mppower(gx, 0.5, 1e-8)
  gx.inv.half = mppower(gx, -0.5, 1e-8)
  A =  mppower(1 / n * gx.half %*% t(xcoef) %*% qmat(n) %*% xcoef %*%
                 gx.half, -0.5, 1e-8)
  M <- A %*% gx.half %*% beta_vector %*% t(beta_vector) %*% gx.half %*% A
  # Obtain eigenvectors and eigenvalues of symmetric M
  gg <- eigen(M, symmetric = T)$vectors
  # Compute directions of central quantile subspace
  vv <- gx.inv.half %*% A %*% gg[, 1:dtau]
  vv2 <- xcoef %*% gx.half %*% A %*% gg[, 1:dtau]

  # Return directions of FCQS
  list(ffun = vv, betax = vv2)
}

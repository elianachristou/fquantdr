#' Scalar-on-function regression
#'
#' \code{sonf} models the relationship between a scalar response variable and
#' one or more functional predictor variables.  The function offers the option
#' to include a second derivative penalty to smooth the estimated regression
#' coefficients.
#'
#' This function models the relationship between a scalar response variable
#' and one or more functional predictor variables using basis functions.  The
#' input `xfd` needs to be a functional object (`fd`).  The function returns
#' the regression coefficients, the functional data objects for these
#' coefficients, the estimated response values, and the design matrix.
#'
#' @param y A numeric vector of length \code{n}, representing the scalar
#'      response.
#' @param xfd A functional data object of class \code{fd}.
#' @param dev2_penalty A logical flag indicating whether to apply a second
#'    derivative penalty (default is `FALSE`).
#' @param lambda A penalty parameter used if `dev2_penalty` is `TRUE`.
#'
#' @return \code{sonf} computes a scalar-on-function regression and returns:
#'    \item{betacoef}{The regression coefficients.}
#'    \item{betafd}{The functional data objects for the beta coefficients.}
#'    \item{yhat}{The \code{n}-dimensional vector of predicted values for `y`.}
#'    \item{X}{The \code{n x (nbasis * p)} design matrix used in the
#'         regression.}
#'
#' @examples
#' # set the parameters
#' n <- 100
#' p <- 5
#' nbasis <- 4
#' nt <- 101
#' time <- seq(0, 1, length.out = nt)
#' eta <- matrix(stats::rnorm(n * p * nbasis), nrow = n, ncol = p * nbasis)
#' # Generate the functional data
#' gen_data <- fundata(n, p, nbasis, time, eta)
#' Xc <- gen_data$xc
#' P <- eigen(stats::cov(eta))$vectors
#' mfpca.scores <- eta %*% P
#' # Prepare the fd object
#' databasis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
#' xcoef_array <- array(0, c(nbasis, n, p))
#' for (k in 1:p) {
#'  xfdk <- fda::smooth.basis(time, t(Xc[, , k]), databasis)$fd
#'  xfdk <- fda::center.fd(xfdk)
#'  xk.coef <- t(xfdk$coef)
#'  xcoef_array[, , k] <- t(xk.coef)
#'  }
#' x.fd <- fda::fd(xcoef_array, databasis)
#' # Generate the model
#' error <- rnorm(n)
#' y <- 3 * mfpca.scores[, 1] + error
#' # Perform scalar-on-function without penalty
#' result_no_penalty <- sonf(y, x.fd, dev2_penalty = FALSE)
#' print("Regression coefficients without penalty:")
#' print(result_no_penalty$betacoef)
#' # Compare y and yhat
#' cor(y, result_no_penalty$yhat)
#' # Perform scalar-on-function regression with penalty
#' lambda <- 0.1
#' result_with_penalty <- sonf(y, x.fd, dev2_penalty = TRUE, lambda = lambda)
#' cor(y, result_with_penalty$yhat)
#'
#' @export
sonf = function(y, xfd, dev2_penalty = FALSE, lambda = NULL) {

  # Check if y is a univariate response
  if (!is.vector(y)) {
    stop(paste("y needs to be a univariate vector."))
  }

  # Check if xfd is a functional object of class 'fd'
  if (!inherits(xfd, "fd")) {
    stop("xfd must be a functional data object of class 'fd'.")
  }

  # Check if xfd is a 3-dimensional array
  if (length(dim(xfd$coef)) != 3) {
    stop("The coefficients of the functional data object 'xfd' must be
         a 3-dimensional array with dimensions (nbasis, n, p), i.e., check
         the dimensions of 'xfd$coef'. If you only have one predictor
         variable, convert xfd$coef into a (nbasis x n x 1) array.")
  }

  # Check if the number of observations for xfd and y agree
  if (length(y) != dim(xfd$coef)[2]) {
    stop("The number of observations of 'y' and the number of observations
         in 'xfd' must agree.")
  }

  # Check if n > p
  if (length(y) <= dim(xfd$coef)[3]) {
    stop(paste("The number of observations of y (", length(y), ") should be
               greater than the number of predictors of xfd
               (", dim(xfd$coef)[3], ").", sep = ""))
  }

  # Check if the number of basis agrees
  if (dim(xfd$coef)[1] != xfd$basis$nbasis) {
    stop(paste("The number of basis functions in the 'xfd$coef' object and
               in the 'xfd$basis' object should agree.  Check
               'dim(xfd$coef)[1]' and 'xfd$basis$nbasis'."))
  }

  # Set the parameters
  n <- dim(xfd$coef)[2]
  nbasis <- dim(xfd$coef)[1]
  p <- dim(xfd$coef)[3]
  if(is.na(p)) p <- 1

  # Extract basis information
  basis <- xfd$basis
  bname <- basis$type

  # Generate the B and DB matrices based on basis type
  if (bname == 'bspline') {
    B <- fda::bsplinepen(basis, 0)  # B_ij = <b_i,b_j>
    DB <- fda::bsplinepen(basis, 2) # DB_ij = <b''_i, b''_j>
  }
  if (bname == 'fourier') {
    B <- diag(nbasis)
    DB <- fda::fourierpen(basis, 2)
  }

  # Prepare the coefficient matrix
  if (p == 1) {
    # Convert to matrix if p = 1
    xcoef <- t(as.matrix(xfd$coef[, , 1]))
  } else {
    # block diagonal of B - p times
    B <- Matrix::bdiag(lapply(as.list(rep("B", p)),
                              function(x) eval(as.symbol(x))))
    # block diagonal of DB - p times
    DB <- Matrix::bdiag(lapply(as.list(rep("DB", p)),
                               function(x) eval(as.symbol(x))))

    xcoef <- aperm(xfd$coef, c(2, 1, 3)) # change dim- n X nbasis X p
    xcoef <- matrix(xcoef, nrow = n)
  }

  # Center the functional predictors and the response
  xcoef_cen <- xcoef - colMeans(xcoef)
  my <- mean(y)
  y_cen <- y - my

  # Construct the design matrix
  X <- as.matrix(xcoef_cen %*% B)

  # Calculate the regression coefficients
  if (dev2_penalty == FALSE) {
    beta_coef <- MASS::ginv(as.matrix(t(X) %*% X)) %*% t(X) %*% y_cen
  } else {
    if (is.null(lambda)) {
      stop("lambda must be provided when dev2_penalty is TRUE")
    }
    beta_coef <- MASS::ginv(as.matrix(t(X) %*% X + lambda * DB)) %*%
      t(X) %*% y_cen
  }

  # Predict the response values
  pred <- X %*% beta_coef + my

  # Construct functional data objects for the beta coefficients
  if (p == 1) {
    betafd <- fda::fd(beta_coef, basis)
  } else {
    betafd <- NULL
    for (j in 1:p) {
      start.ind <- (j - 1) * nbasis + 1
      end.ind <- j * nbasis
      betafd[[j]] <- fda::fd(beta_coef[start.ind:end.ind], basis)
    }
  }
  return(list(betacoef = beta_coef, betafd = betafd, yhat = pred, X = X))
}

#' Scalar-on-function regression
#'
#' \code{sonf} models the relationship between a scalar response variable and
#' one or more functional predictor variables.  The function offers the option
#' to include a second derivative penalty to smooth the estimated regression
#' coefficients.
#'
#' This function models the relationship between a scalar response variable
#' and one or more functional predictor variables using basis functions.  The
#' input `xfd' needs to be a functional object ('fd').  The function returns
#' the regression coefficients, the functional data objects for these
#' coefficients, the estimated response values, and the design matrix.
#'
#' @param y A numeric vector of length `n', representing the scalar response.
#' @param xfd A functional data object of class `fd'.
#' @param dev2_penalty A logical flag indicating whether to apply a second
#'    derivative penalty (default is `FALSE').
#' @param lambda A penalty parameter used if `dev2_penalty' is `TRUE'.
#'
#' @return \code{sonf} computes a scalar-on-function regression and returns:
#'    \item{betacoef}{The regression coefficients.}
#'    \item{betafd}{The functional data objects for the beta coefficients.}
#'    \item{yhat}{The \code{n}-dimensional vector of predicted values for `y'.}
#'    \item{X}{The \code{n x (nbasis * p)} design matrix used in the
#'         regression.}
#'
#' @examples
#' # Example 1
#' # Set the parameters
#' n <- 100
#' nbasis <- 10
#' p <- 3
#' # Create a B-spline bsis
#' basis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis,
#'                norder = nbasis)
#' # Generate random coefficients for the functional data object
#' coef_matrix <- array(rnorm(n * nbasis * p), dim = c(nbasis, n, p))
#' fdobj <- fda::fd(coef_matrix, basis)
#' # Create a scalar response vector
#' y <- rnorm(n)
#' # Perform scalar-on-function regression without penalty
#' result_no_penalty <- sonf(y, fdobj, dev2_penalty = FALSE)
#' print("Regression coefficients without penalty:")
#' print(result_no_penalty$betacoef)
#' # Perform scalar-on-function regression with penalty
#' lambda <- 0.1
#' result_with_penalty <- sonf(y, fdobj, dev2_penalty = TRUE,
#'           lambda = lambda)
#' print("Regression coefficients with penalty:")
#' print(result_with_penalty$betacoef)
#' print("Predicted values:")
#' print(result_with_penalty$yhat)
#'
#' # Example 2
#' # set the parameters
#' n <- 100
#' p <- 5
#' nbasis <- 4
#' nt <- 101
#' time <- seq(0, 1, length.out = nt)
#' eta <- matrix(stats::rnorm(n * p * nbasis), nrow = n,
#'     ncol = p * nbasis)
#' # Generate the functional data
#' gen_data <- fundata(n, p, nbasis, time, eta)
#' Xc <- gen_data$xc
#' P <- eigen(stats::cov(eta))$vectors
#' mfpca.scores <- eta %*% P
#' # Prepare the fd object
#' databasis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis,
#' norder = nbasis)
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
#'
#' @export
sonf = function(y, xfd, dev2_penalty = FALSE, lambda = NULL) {

  y <- as.matrix(y) # double check!

  # Check if y is univariate response
  if (dim(y)[2] > 1) {
    stop(paste("y needs to be a univariate response. y is a",
               dim(y)[2], "-dimensional response in this case."))
  }

  # Check if xfd is a functional object of class 'fd'
  if (!inherits(xfd, "fd")) {
    stop("xfd must be a functional data object of class 'fd'")
  }

  # Check if xfd is a 3-dimensional array
  if (length(dim(xfd$coef)) != 3) {
    stop("xfd$coef must be an array with dimensions (nbasis, n, p).")
  }

  # Check if the number of observations for xfd and y agree
  if (length(y) != dim(xfd$coef)[2]) {
    stop("The length of y and the number of observations in xfd must agree.")
  }

  # Check if n > p
  if (length(y) <= dim(xfd$coef)[3]) {
    stop(paste("number of observations of y (", length(y), ") should be
               greater than the number of predictors of xfd
               (", dim(xfd$coef)[3], ").", sep = ""))
  }

  # Extract basis information
  basis <- xfd$basis
  nbasis <- basis$nbasis
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

  # Set the parameters
  tmp <- dim(xfd$coef)
  n <- tmp[2]
  nbasis <- tmp[1]
  p <- tmp[3]
  if(is.na(p)) p <- 1

  # Check the number of functional predictors
  if (p != floor(p) | p <= 0) {
    stop("The number of functional predictors (p) must be a positive
         integer.")
  }

  # Prepare the coefficient matrix
  if (p == 1) {
    # Convert to. matrix if p = 1
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
    #(xvec[1,] == as.vector(xcoef[1,,]))
  }

  # Center the functional predictors and the response
  xcoef_cen <- xcoef-colMeans(xcoef)
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

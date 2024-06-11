#' Scalar-on-function regression
#'
#' \code{sonf} estimates the relationship between a scalar response and
#' functional predictors using basis functions, with options for applying
#' a second derivative penalty or not.
#'
#' This function performs scalar-on-function regression, predicting a scalar
#' response \code{y} using one or more functional predictors \code{xfd}. It
#' centralizes the data, constructs a design matrix based on the different type
#' of basis functions, and calculates the regression coefficients with or without
#' a second derivative penalty. The function returns the regression coefficients,
#' functional data objects for these coefficients, estimated response value,
#' and the designed matrix X.
#'
#' @param y A number vector of length n, representing the scalar response variable
#' @param xfd A functional data object
#' @param dev2_penalty A logical flag indicating whether to apply a second
#'    derivative penalty.
#' @param lambda A penalty parameter (used if \code{dev2_penalty} is 'TRUE')
#'
#' @return \code{sonf} computes a scalar-on-function regression and returns:
#'    \itemize{
#'    \item{beta_coef: }{The regression coefficients.}
#'
#'    \item{beta_fd: }{The functional data objects for the beta coefficients.}
#'
#'    \item{yhat: }{The predicted values of \code{y}}
#'
#'    \item{X: }{The designed matrix used in the regression}}
#'
#' @examples
#' # example code
#'
#' @export
sonf = function(y, xfd, dev2_penalty = FALSE, lambda=NULL) {

  y <- as.matrix(y) # double check!

  # compatibility checks
  # checks if y is univariate
  if (dim(y)[2] > 1) {
    stop(paste("y needs to be a univariate response. y is a", dim(y)[2], "-dimensional response in this case."))
  }
  # checks if xfd is multivariate
  if (!inherits(xfd, "fd")) {
    stop("xfd must be a functional data object of class 'fd'")
  }
  # checks if the coefficients of xfd is an array
  if (length(dim(xfd$coef)) != 3) {
    stop("xfd$coef must be an array with dimensions (nbasis, n, p)")
  }
  # checks if the number of observations for xfd and y agree
  if (length(y) != dim(xfd$coef)[2]) {
    stop("The length of y and the number of observations in xfd must agree")
  }
  # checks if n > p
  if (length(y) <= dim(xfd$coef)[3]) {
    stop(paste("number of observations of y (", length(y), ") should be greater than the number of predictors of xfd (", dim(xfd$coef)[3], ").", sep = ""))
  }

  basis <- xfd$basis
  nbasis <- basis$nbasis
  bname <- basis$type
  if (bname == 'bspline') {
    B <- bsplinepen(basis, 0)  # B_ij = <b_i,b_j>
    DB <- bsplinepen(basis, 2) # DB_ij = <b''_i, b''_j>
  }
  if (bname == 'fourier') {
    B <- diag(nbasis)
    DB <- fourierpen(basis, 2)
  }
  tmp <- dim(xfd$coef)
  n <- tmp[2]
  nbasis <- tmp[1]
  p <- tmp[3]
  if(is.na(p)) p <- 1

  if (!is.numeric(p) || p <= 0) {
    stop("The number of functional predictors (p) must be a positive integer")
  }

  if (p == 1) {
    xcoef <- t(as.matrix(xfd$coef[, , 1])) # because originally it is an array
  } else {
    # block diagonal of B - p times
    B <- Matrix::bdiag(lapply(as.list(rep("B", p)), function(x) eval(as.symbol(x))))
    # block diagonal of DB - p times
    DB <- Matrix::bdiag(lapply(as.list(rep("DB", p)), function(x) eval(as.symbol(x))))

    xcoef <- aperm(xfd$coef, c(2, 1, 3)) # change dim- n X nbasis X p
    xcoef <- matrix(xcoef, nrow = n)
    #(xvec[1,] == as.vector(xcoef[1,,]))
  }
  xcoef_cen <- xcoef-colMeans(xcoef)
  my <- mean(y)
  y_cen <- y - my
  X <- as.matrix(xcoef_cen %*% B)
  if (dev2_penalty == FALSE) {
    beta_coef <- ginv(as.matrix(t(X) %*% X)) %*% t(X) %*% y_cen
  } else {
    if (is.null(lambda)) {
      stop("lambda must be provided when dev2_penalty is TRUE")
    }
    beta_coef <- ginv(as.matrix(t(X) %*% X + lambda * DB)) %*% t(X) %*% y_cen
  }
  pred <- X %*% beta_coef + my
  if (p == 1) {
    betafd <- fd(beta_coef, basis)
  } else {
    betafd <- NULL
    for (j in 1:p) {
      start.ind <- (j - 1) * nbasis + 1
      end.ind <- j * nbasis
      betafd[[j]] <- fd(beta_coef[start.ind:end.ind], basis)
    }
  }
  return(list(beta_coef = beta_coef, betafd = betafd, yhat = pred, X = X))
}

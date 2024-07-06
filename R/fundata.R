#' Functional data
#'
#' \code{fundata} generates functional data based on Fourier basis functions.
#'
#' This function generates functional data for a given number of observations
#' ('n'), functional predictors ('p'), basis functions ('nbasis'), and time
#' points ('time').  It utilizes a specified coefficient matrix ('eta') to
#' create the functional predictors.  Using these basis functions, it computes
#' the original functional predictors ('x') by multiplying the coefficients
#' with the basis functions. The function then centers these functional
#' predictors ('xc') by subtracting the mean across observations for each time
#' point.  The output is a list containing both the original and centered
#' functional predictors.
#'
#' @param n The number of observations, i.e., sample size.
#' @param p The number of functional predictors.
#' @param nbasis The number of Fourier basis functions.
#' @param time A vector of time points at which the functional data is
#'     evaluated.
#' @param eta A matrix of coefficients with dimensions \code{n} by
#'     \code{p * nbasis}.
#'
#' @return \code{fundata} generates functional data based on Fourier
#'    basis functions and returns:
#'    \itemize{
#'        \item \code{x}: The original functional predictors, a
#'            \code{n * nt * p} array, where \code{nt} denotes the number
#'            of time points.
#'        \item \code{xc}: The centered functional predictors, a
#'            \code{n * nt * p} array, where \code{nt} denotes the number
#'            of time points.
#'    }
#'
#' @examples
#' # Example 1
#' # set the parameters
#' n <- 100
#' p <- 5
#' nbasis <- 3
#' time <- seq(0, 1, length.out = 101)
#' eta <- matrix(stats::rnorm(n * p * nbasis), nrow = n, ncol = p * nbasis)
#' # create the functional predictors
#' result <- fundata(n, p, nbasis, time, eta)
#' x <- result$x
#' xc <- result$xc
#' # plot the first functional predictor for illustration
#'  fda::matplot(time, t(x[, , 1]), type = "l", lty = 1, col = 1:n,
#'  xlab = "Time", ylab = "Value", main = paste("Functional Predictor", 1))
#'
#' # Example 2
#' # set the parameters
#' n <- 100
#' p <- 5
#' nbasis <- 4
#' time <- seq(0, 1, length.out = 101)
#' SigmaCov <- matrix(0, nrow = p * nbasis, ncol = p * nbasis)
#' for (j in 1:p) {
#'  index.j <-(((j - 1) * nbasis + 1):(j * nbasis))
#'  diag(SigmaCov[index.j, index.j]) <- c(2, 1, 1/2, 1/4)
#' }
#' eta <- mvtnorm::rmvnorm(n, mean = rep(0, p * nbasis), sigma = SigmaCov)
#' # create the functional predictors
#' result <- fundata(n, p, nbasis, time, eta)
#' x <- result$x
#' xc <- result$xc
#' # plot the first functional predictor for illustration
#' fda::matplot(time, t(x[, , 1]), type = "l", lty = 1, col = 1:n,
#'     xlab = "Time", ylab = "Value", main = paste("Functional Predictor", 1))
#'
#' @export
fundata <- function(n, p, nbasis, time, eta) {

  # Check if n is a single integer number
  if (length(n) != 1 | n != round(n)) {
    stop("Parameter 'n' must be a single integer number.")
  }

  # Check if n is positive and greater than p
  if (n <= 0 | n <= p) {
    stop("Parameter 'n' must be a positive integer and greater than 'p'.")
  }

  # Check that p is a single number
  if (length(p) != 1) {
    stop("Parameter 'p' must be a single number.")
  }

  # Check that p is a positive integer
  if (p != round(p) | p <=0) {
    stop("Parameter 'p' must be a positive integer number.")
  }

  # Check that nbasis is a single number
  if (length(nbasis) != 1) {
    stop("Parameter 'nbasis' must be a single number.")
  }

  # Check that nbasis is a positive integer
  if (nbasis != round(nbasis) | nbasis <=0) {
    stop("Parameter 'nbasis' must be a positive integer number.")
  }

  # Check if time is a vector
  if (!is.vector(time) | length(time) == 1) {
    stop("time is a vector of length more than 1.")
  }

  # Check if the dimension of eta is n * (p * nbasis)
  if (!is.matrix(eta) | nrow(eta) != n | ncol(eta) != (p * nbasis)) {
    stop("Parameter 'eta' must be a numeric matrix with dimensions
         n x (p * nbasis).")
  }

  # define the parameters
  nt <- length(time)
  x <- array(0, dim = c(n, nt, p))
  xc <- array(0, dim = c(n, nt, p))

  ## Create a Fourier basis with nbasis functions over the interval [0, 1]
  f.ans <- fda::create.fourier.basis(rangeval = c(0, 1), nbasis = nbasis,
                                     dropind = 1)
  st <- as.matrix(fda::eval.basis(time, f.ans))

  # when nbasis is even, the dimension of st is n x nbasis, but when
  # nbasis is odd, the dimension of st is n x (nbasis - 1)
  nbasis <- dim(st)[2]

  ## Generate functional predictors
  # x and xc are the original and centered functional predictors
  for(j in 1:p) {
    index.j <- (((j - 1) * nbasis + 1):(j * nbasis))
    x[, , j] <- eta[, index.j] %*% t(st)
    xc[, , j] <- x[, , j] - matrix(rep(apply(x[, , j], 2, mean), n),
                                   nrow = n, byrow = T)
  }
  return(list(x = x, xc = xc))
}

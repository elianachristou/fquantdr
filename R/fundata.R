#' Functional data
#'
#' \code{fundata} generates functional data based on Fourier basis functions.
#'
#' This function generates functional data for a given number of observations
#' ('n'), functional predictors ('p'), basis functions ('nbasis'), and times
#' points ('time').  It utilizes a specified coefficient matrix ('eta') to
#' create the functional predictors.  Using these basis functions, it computes
#' the original functional predictors ('g') by multiplying the coefficients
#' with the basis functions. The function then centers these functional
#' predictors ('cg') by subtracting the mean across observations for each time
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
#'        \item \code{g}: The original functional predictors, a
#'            \code{n * nt * p} array, where \code{nt} denotes the number
#'            of time points.
#'        \item \code{cg}: The centered functional predictors, a
#'            \code{n * nt * p} array, where \code{nt} denotes the number
#'            of time points.
#'    }
#'
#' @examples
#' # Example 1
#' n <- 100
#' p <- 5
#' nbasis <- 3
#' time <- seq(0, 1, length.out = 101)
#' eta <- matrix(stats::rnorm(n * p * nbasis), nrow = n, ncol = p * nbasis)
#' result <- fundata(n, p, nbasis, time, eta)
#' # original and centered functional predictors
#' g <- result$g
#' cg <- result$cg
#' # plot the first functional predictor for illustration purposes
#'  fda::matplot(time, t(g[, , 1]), type = "l", lty = 1, col = 1:n,
#'  xlab = "Time", ylab = "Value", main = paste("Functional Predictor", 1))
#'
#' # Example 2
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
#' result <- fundata(n, p, nbasis, time, eta)
#' # plot the first functional predictor for illustration purposes
#' fda::matplot(time, t(result$g[, , 1]), type = "l", lty = 1, col = 1:n,
#'     xlab = "Time", ylab = "Value", main = paste("Functional Predictor", 1))
#'
#' @export
fundata <- function(n, p, nbasis, time, eta) {

  # compatibility checks
  # checks if n is an integer and n > p
  if (length(n) != 1 | n != round(n)) {
    stop("Parameter 'n' must be a single integer number.")
  }

  if (n <= 0 | n <= p) {
    stop("Parameter 'n' must be a positive integer and greater than 'p'.")
  }

  # checks that p is a single number
  if (length(p) != 1) {
    stop("Parameter 'p' must be a single number.")
  }

  # checks that p is a positive integer
  if (p != round(p) | p <=0) {
    stop("Parameter 'p' must be integer and positive number.")
  }

  # checks that nbasis is a single number
  if (length(nbasis) != 1) {
    stop("Parameter 'nbasis' must be a single number.")
  }

  # checks that nbasis is a positive integer
  if (nbasis != round(nbasis) | nbasis <=0) {
    stop("Parameter 'nbasis' must be integer and positive number.")
  }

  # check if time is a vector
  if (!is.vector(time) | length(time) == 1) {
    stop("time is a vector of length more than 1.")
  }

  # checks if the dimension of eta is n * (p * nbasis)
  if (!is.matrix(eta) | nrow(eta) != n | ncol(eta) != (p * nbasis)) {
    stop("Parameter 'eta' must be a numeric matrix with dimensions
         n x (p * nbasis).")
  }

  # define the parameters
  nt <- length(time)
  g <- array(0, dim = c(n, nt, p))
  cg <- array(0, dim = c(n, nt, p))

  ## Create a Fourier basis with nbasis functions over the interval [0, 1]
  f.ans <- fda::create.fourier.basis(rangeval = c(0, 1), nbasis = nbasis,
                                     dropind = 1)
  st <- as.matrix(fda::eval.basis(time, f.ans))

  # when nbasis is even, the dimension of st is n x nbasis, but when
  # nbasis is odd, the dimension of st is n x (nbasis - 1)
  nbasis <- dim(st)[2]

  ## Generate functional predictors
  # g and cg are the original and centered functional predictors
  for(j in 1:p) {
    index.j <- (((j - 1) * nbasis + 1):(j * nbasis))
    g[, , j] <- eta[, index.j] %*% t(st)
    cg[, , j] <- g[, , j] - matrix(rep(apply(g[, , j], 2, mean), n),
                                   nrow = n, byrow = T)
  }
  return(list(g = g, cg = cg))
}

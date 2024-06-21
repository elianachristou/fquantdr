#' Functional data
#'
#' \code{fundata} generates functional data based on Fourier basis functions
#'
#' This function generates functional data for a given number of observations
#' ('n'), functional predictors ('p'), basis functions ('q'), and times points
#' ('t').  It utilizes a specified coefficient matrix ('eta') to create the
#' functional predictors.  Using these basis functions, it computes the
#' original functional predictors ('g') by multiplying the coefficients
#' with the basis functions. The function then centers these functional
#' predictors ('cg') by subtracting the mean across observations for each
#' time point.  The output is a list containing both the original and centered
#' functional predictors.
#'
#' @param n The number of observations, i.e., sample size.
#' @param p The number of functional predictors.
#' @param q The number of Fourier basis functions.
#' @param t A vector of time points at which the functional data is evaluated.
#' @param eta A matrix of coefficients with dimensions \code{n} by \code{p * q}.
#'
#' @return \code{fundata} generates functional data based on Fourier
#'    basis functions and returns:
#'    \itemize{
#'        \item \code{g}: The original functional predictors, an
#'            \code{n * length(t) * p} array.
#'        \item \code{cg}: The centered functional predictors, an
#'            \code{n * length(t) * p} array.
#'    }
#'
#' @examples
#' # Example 1
#' n <- 100
#' p <- 5
#' q <- 4
#' t <- seq(0, 1, length.out = 101)
#' eta <- matrix(stats::rnorm(n * p * q), nrow = n, ncol = p * q)
#' result <- fundata(n, p, q, t, eta)
#' # original functional predictors
#' g <- result$g
#' # centered functional predictors
#' cg <- result$cg
#' # plot the first functional predictor for illustration purposes
#'  fda::matplot(t, t(g[, , 1]), type = "l", lty = 1, col = 1:n, xlab = "Time",
#'  ylab = "Value", main = paste("Functional Predictor", 1))
#'
#' # Example 2
#' n <- 100
#' p <- 5
#' q <- 4
#' t <- seq(0, 1, length.out = 101)
#' SigmaCov <- matrix(0, nrow = p * q, ncol = p * q)
#' for (j in 1:p) {
#'  index.j <-(((j - 1) * q + 1):(j * q))
#'  diag(SigmaCov[index.j, index.j]) <- c(2, 1, 1/2, 1/4)
#' }
#' eta <- mvtnorm::rmvnorm(n, mean = rep(0, p * q), sigma = SigmaCov)
#' result <- fundata(n, p, q, t, eta)
#' # plot the first functional predictor for illustration purposes
#' fda::matplot(t, t(result$g[, , 1]), type = "l", lty = 1, col = 1:n,
#'     xlab = "Time", ylab = "Value", main = paste("Functional Predictor", 1))
#'
#' @export
fundata <- function(n, p, q, t, eta) {

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

  # checks that q is a single number
  if (length(q) != 1) {
    stop("Parameter 'q' must be a single number.")
  }

  # checks that p is a positive integer
  if (p != round(p) | p <=0) {
    stop("Parameter 'p' must be integer and positive number.")
  }

  # checks that q is a positive integer
  if (q != round(q) | q <=0) {
    stop("Parameter 'q' must be integer and positive number.")
  }

  # check if t is a vector
  if (is.vector(t) | length(t) == 1) {
    stop("t is a vector of length more than 1.")
  }

  # checks if the dimension of eta is n * (p * q)
  if (!is.matrix(eta) | nrow(eta) != n | ncol(eta) != (p * q)) {
    stop("Parameter 'eta' must be a numeric matrix with dimensions n x (p * q).")
  }

  # define the parameters
  Time <- length(t)
  g <- array(0, dim = c(n, Time, p))
  cg <- array(0, dim = c(n, Time, p))

  ## Create a Fourier basis with q basis functions over the interval [0, 1]
  f.ans <- fda::create.fourier.basis(rangeval = c(0, 1), nbasis = q,
                                     dropind = 1)
  st <- as.matrix(fda::eval.basis(t, f.ans))

  # when q is even, the dimension of st is n x q, but when q is odd,
      # the dimension of st is n x (q - 1)
  q <- dim(st)[2]

  ## Generate functional predictors
  # g and cg are the original and centered functional predictors
  for(j in 1:p) {
    index.j <- (((j - 1) * q + 1):(j * q))
    g[, , j] <- eta[, index.j] %*% t(st)
    cg[, , j] <- g[, , j] - matrix(rep(apply(g[, , j], 2, mean), n), nrow = n, byrow = T)
  }
  return(list(g = g, cg = cg))
}

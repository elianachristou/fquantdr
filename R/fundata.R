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
#'        \item \code{g}: the original functional predictors
#'        \item \code{cg}: the centered functional predictors
#'    }
#'
#' @references Wang, G., Liu, S., Han, F., and Di, C.-Z. (2022). Robust functional
#'    principal component analysis via a functional pairwise spatial sign operator.
#'    \emph{Biometrics} https://doi.org/10.1111/biom.13695.
#'
#' @import stats
#' @examples
#' # Example 1
#' n <- 100
#' p <- 5
#' q <- 4 #(double check!)
#' t <- seq(0, 1, length.out = 101)
#' eta <- matrix(rnorm(n * p * q), n, p * q)
#'
#' result <- fundata(n, p, q, t, eta)
#' g <- result$g # Original functional predictors
#' cg <- result$cg # Centered functional predictors
#'
#' @export
fundata <- function(n, p, q, t, eta) {
  # define the parameters
  Time <- length(t)
  g <- array(0, dim = c(n, Time, p))
  cg <- array(0, dim = c(n, Time, p))

  # compatibility checks
  # checks if the dimension of eta is n * (p * q)
  if (!is.matrix(eta) || nrow(eta) != n || ncol(eta) != (p * q)) {
    stop("Parameter 'eta' must be a numeric matrix with dimensions n x (p * q).")
  }

  # checks if n is an interger and n > p
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n <= p) {
    stop("Parameter 'n' must be a positive integer and greater than 'p'.")
  }

  ## Create a Fourier basis with q basis functions over the interval [0, 1]
  f.ans <- fda::create.fourier.basis(rangeval = c(0, 1), nbasis = q,
                                     dropind = 1)
  st <- as.matrix(fda::eval.basis(t, f.ans))

  ## Generate functional predictors
  # g and cg are the original and centered functional predictors
  for(j in 1:p) {
    index.j <- (((j - 1) * q + 1):(j * q))
    g[, , j] <- eta[, index.j] %*% t(st)
    cg[, , j] <- g[, , j] - matrix(rep(apply(g[, , j], 2, mean), n), nrow = n, byrow = T)
  }
  return(list(g = g, cg = cg))
}

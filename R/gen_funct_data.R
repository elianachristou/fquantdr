#' function that generates noiseless functional data using Fourier Basis functions
#'
#' \code{gen_funct_data} generates noiseless functional data using Fourier
#' basis functions
#'
#' This function initializes arrays for the original and centered functional
#' predictors and construct Fourier basis functions over the interval [0, 1].
#' For each predictor, it computes the functional data by multiplying the
#' relevant scores with the evaluated basis functions and then centers this data
#' by subtracting the mean of each time point across all observation.
#'
#' @param n Sample size, i.e., the number of observations
#' @param p the number of predictors
#' @param q the number of Fourier basis functions used to generate the
#'    functional data
#' @param t time points
#' @param eta A matrix of scores of size pq * pq. These scores can be Gaussian
#'    or elliptical
#'
#' @return \code{gen_funct_data} generates the functional data using Fourier
#'    basis functions and returns:
#'    \itemize{
#'    \item{g: }{the original functional predictors}
#'    \item{cg: }{the centered functional predictors}
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
#' p <- 3
#' q <- 5
#' t <- seq(0, 1, length.out = 101)
#' eta <- matrix(rnorm(n * p * q), n, p * q)
#'
#' result <- gen_funct_data(n, p, q, t, eta)
#' g <- result$g # Original functional predictors
#' cg <- result$cg # Centered functional predictors
#'
#' @export
gen_funct_data <- function(n, p, q, t, eta) {
  Time <- length(t)
  g <- array(0, dim = c(n, Time, p))
  cg <- array(0, dim = c(n, Time, p))

  ## Fourier basis
  f.ans <- fda::create.fourier.basis(rangeval = c(0, 1), nbasis = q, dropind = 1)
  st <- as.matrix(fda::eval.basis(t, f.ans))

  ## g and cg are original and centered functional predictors
  for(j in 1:p) {
    index.j <- (((j - 1) * q + 1):(j * q))
    g[, , j] <- eta[, index.j] %*% t(st)
    cg[, , j] <- g[, , j] - matrix(rep(apply(g[, , j], 2, mean), n), nrow = n, byrow = T) }
  return(list(g = g, cg = cg))
}

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
#'





gen_funct_data <- function(n, p, q, t, eta) {
  Time <- length(t)
  g <- array(0, dim = c(n, Time, p))
  cg <- array(0, dim = c(n, Time, p))

  ## Fourier basis
  f.ans <- create.fourier.basis(rangeval = c(0, 1), nbasis = q, dropind = 1)
  st <- as.matrix(eval.basis(t, f.ans))

  ## g and cg are original and centered functional predictors
  for(j in 1:p) {
    index.j <- (((j - 1) * q + 1):(j * q))
    g[, , j] <- eta[, index.j] %*% t(st)
    cg[, , j] <- g[, , j] - matrix(rep(apply(g[, , j], 2, mean), n), nrow = n, byrow = T) }
  return(list(g = g, cg = cg))
}

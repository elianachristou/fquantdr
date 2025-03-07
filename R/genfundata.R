#' Functional data
#'
#' \code{genfundata} generates functional data using either B-spline or Fourier
#' basis functions.
#'
#' This function constructs functional data using basis expansions.  The
#' coefficients for the basis functions are randomly generated, and an
#' orthogonalization step is applied using a Q-matrix.  Functional principal
#' component analysis (FPCA) is then performed to obtain a low-dimensional
#' representation of the data.
#'
#' @param n The number of observations, i.e., sample size.
#' @param p The number of functional predictors.
#' @param nbasis The number of basis functions.
#' @param tt A vector of time points at which the functional data is
#'     evaluated.
#' @param basisname A character string for the type of basis functions to use,
#'     either 'bspline' or 'fourier'.  Default is 'bspline'.
#'
#' @return \code{genfundata} returns a list containing:
#'    \itemize{
#'        \item \code{X}: An array of generated functional data of dimensions
#'            \code{n * nt * p} array, where \code{nt} denotes the number
#'            of time points.
#'        \item \code{mfpca.scores}: Matrix of FPCA scores.
#'    }
#'
#' @examples
#' n <- 50
#' p <- 2
#' nbasis <- 5
#' tt <- seq(0, 1, length.out = 100)
#' data <- genfundata(n, p, nbasis, tt, 'bspline')
#' str(data)
#'
#' @export


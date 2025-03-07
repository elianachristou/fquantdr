#' Functional Principal Component Analysis
#'
#' \code{fpca} performs Functional Principal Component Analysis (FPCA) on a given
#' functional object.
#'
#' This function performs FPCA on a given function object and returns the
#' eigenvalues and eigenfunctions of the covariance matrix.
#'
#' @param ftn A functional data object containing coefficients and basis information.
#' @param basisname A character string specifying the type of basis function to use.
#'     Options are 'bspline' or 'fourier'.
#'
#' @return \code{fpca} returns a list containing:
#'    \itemize{
#'        \item \code{pred}: Principal component scores.
#'        \item \code{eval}: Eigenvalues of the covariance operator.
#'        \item \code{mat}: Transformation matrix for principal components.
#'    }
#'
#' @examples
#'
#' @export
fpca <- function(ftn, basisname) {
  # Extract coefficient array from input 'ftn' and determine its dimensions
  temp <- ftn$coef
  n <- dim(temp)[2]
  p <- dim(temp)[3]
  nt <- dim(temp)[1]

  if(is.na(p)) p <- 1
  xcoef <- temp

  basis <- ftn$basis
  # Generate basis matrix based on specified basis type
  if(basisname == 'bspline') {
    GB <- bsplinepen(basis, Lfdobj = 0)
  }
  else if(basisname == 'fourier') {
    GB <- fourierpen(basis, Lfdobj = 0)
  }

  # Centering matrix for FPCA
  one <- matrix(1, n, 1)
  Q <- diag(n) - one %*% t(one) / n
  # Compute the square root of the penalty matrix
  B.half <- matpower(GB, 0.5)

  if(p==1) {
    Sigma <- B.half %*% xcoef %*% Q %*% t(xcoef) %*% B.half / n
    egn <- eigen(Sigma, sym = TRUE)
    B.inv.half <- matpower(GB, -0.5)
    pred <- Q %*% t(xcoef) %*% GB %*% B.inv.half %*% egn$vec
    out <- list(pred = pred, eval = egn$val, mat = B.inv.half %*% egn$vec)
  } else if(p > 1) { # BX is the same for now...
    M.half <- B.half %*% xcoef[, , 1] %*% Q
    B.inv.half <- matpower(GB, -0.5)
    D.inv.half <- B.inv.half
    BB <- GB %*% xcoef[, , 1] %*% Q
    for(j in 2:p) {
      D.inv.half <- as.matrix(bdiag2(D.inv.half, B.inv.half))
      M.half <- rbind(M.half, B.half %*% xcoef[, , j] %*% Q)
      BB <- rbind(BB, GB %*% xcoef[, , j] %*% Q)
    }
    DD <- M.half
    M.half <- M.half %*% t(M.half) / n
    egn <- eigen(M.half, sym = TRUE)
    pred <- t(BB) %*% D.inv.half %*% egn$vec
    out <- list(pred = pred, eval = egn$val, mat = D.inv.half %*% egn$vec)
  }
  return(out)
}

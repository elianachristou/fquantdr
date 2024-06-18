#' mfsir: FSIR that can incorporate multivariate functional predictor

#' \code{mfsir} description

#' This code (description)

#' (X can be p-dimensional and must be centered)

#' @param Xc is a (n x t x p) functional element - array
#' @param y is the response
#' @param H is the number of slices
#' @param gx is the gram matrix of the basis functions
#' @param nbasis is the number of basis functions

#' @return the eigenvectors and the sufficient predictors

#' @noRd

#' @examples
#'


mfsir <- function(Xc, y, H, nbasis) {
  n <- dim(Xc)[1]
  p <- dim(Xc)[3]

  databasis <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
  xfd.coef <- numeric() #retrieves the coodinates as an n x (p * q) matrix
  for (k in 1:p) {
    xfdk <- fda::smooth.basis(t, t(Xc[, , k]), databasis)$fd
    xfdk <- fda::center.fd(xfdk)
    xk.coef <- t(xfdk$coef)
    xfd.coef <- cbind(xfd.coef, xk.coef)
  }

  #compute block diagonal gram matrix gx
  gi <- gramatrix(nbasis, databasis)

  gx <- matrix(0, nrow = p * (nbasis), ncol = p * (nbasis))
  for (i in 1:p) {
    index <- ((i - 1) * (nbasis) + 1):(i * (nbasis))
    gx[index, index] <- gi
  }

  sigmaxx <- (t(xfd.coef) %*% xfd.coef) / n
  gx.sqrt <- mppower(gx, 1/2, 10e-7)
  store <- gx.sqrt %*% sigmaxx %*% gx.sqrt
  store.inv.sqrt <- rigpower(store, -1/2, 0.05)
  mucoef <- apply(xfd.coef, 2, mean)
  ydis <- discretize (y, 1:H)
  prob <- slprob(ydis, 1:H)
  avg.slav <- slav(xfd.coef, ydis, 1:H)
  Lambda1 <- matrix(0, p * nbasis, p * nbasis)
  for(i in 1:H) {
    Lambda1 <- Lambda1 + prob[i] * avg.slav[i, ] %*% t(avg.slav[i, ])}
  sigmaxxy <- Lambda1 - mucoef %*% (t(mucoef))
  mat <- store.inv.sqrt %*% gx.sqrt
  out <- mat %*% sigmaxxy %*% t(mat)
  out <- symmetry(out)
  eigendecom <- eigen(out)
  eigvalues <- eigendecom$values
  phi <- eigendecom$vec
  gx.inv.sqrt <- rigpower(gx, -1/2, 0.05)
  betas <- gx.inv.sqrt %*% store.inv.sqrt %*% phi
  return(list(phi = phi, betas = betas, eigvalues = eigvalues, xfd.coef = xfd.coef, gx = gx))
}

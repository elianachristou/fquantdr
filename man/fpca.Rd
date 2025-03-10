% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fpca.R
\name{fpca}
\alias{fpca}
\title{Functional Principal Component Analysis}
\usage{
fpca(ftn, basisname)
}
\arguments{
\item{ftn}{A functional data object containing coefficients and basis information.
The coefficients is a 3D array (q x n x p), where \code{q} is the number of basis,
\code{n} is the sample size, and \code{p} is the number of functional variables.}

\item{basisname}{A character string specifying the type of basis function to use.
Options are 'bspline' or 'fourier'.}
}
\value{
\code{fpca} returns a list containing:
\itemize{
\item \code{pred}: Principal component scores.
\item \code{eval}: Eigenvalues of the covariance operator.
\item \code{mat}: Transformation matrix for principal components.
}
}
\description{
\code{fpca} computes FPCA by first extracting the coefficient matrix from the
functional data object, then constructing the appropriate penalty matrix based
on the specified basis type.  The method applied eigen decomposition to obtain
functional principal components.
}
\details{
This function performs FPCA on a given function object and returns the
eigenvalues and eigenfunctions of the covariance matrix.
}
\examples{
# set the parameters
n <- 100
p <- 5
nbasis <- 4
nt <- 101
time <- seq(0, 1, length.out = nt)
eta <- matrix(stats::rnorm(n * p * nbasis), nrow = n, ncol = p * nbasis)
# Generate the functional data
gen_data <- fundata(n, p, nbasis, time, eta)
Xc <- gen_data$xc
P <- eigen(stats::cov(eta))$vectors
mfpca.scores <- eta \%*\% P
# Prepare the fd object
databasis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
xcoef_array <- array(0, c(nbasis, n, p))
for (k in 1:p) {
 xfdk <- fda::smooth.basis(time, t(Xc[, , k]), databasis)$fd
 xfdk <- fda::center.fd(xfdk)
 xk.coef <- t(xfdk$coef)
 xcoef_array[, , k] <- t(xk.coef)
 }
x.fd <- fda::fd(xcoef_array, databasis)

}

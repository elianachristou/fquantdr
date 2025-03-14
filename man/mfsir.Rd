% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mfsir.R
\name{mfsir}
\alias{mfsir}
\title{Functional Sliced Inverse Regression}
\usage{
mfsir(X, y, H, nbasis = 4, norder = 4)
}
\arguments{
\item{X}{A 3-dimensional array (\code{n x nt x p}), where n is the number
of observations, nt is the number of time points, and p is the number
of predictor variables.}

\item{y}{A numeric vector of length \code{n} representing the scalar response
variable.}

\item{H}{An integer specifying the number of slices for the response variable;
must be at least 2 and less than \code{n}.}

\item{nbasis}{An integer specifying the number of basis functions for functional
smoothing.  Currently, the only option is to use B-spline basis. Default
value is 4.}

\item{norder}{An integer specifying the order of B-splines, which is one
higher than their degree.  The default of 4 gives cubic splines.}
}
\value{
\code{mfsir} computes the new sufficient predictors and returns
\item{xcoef}{A \code{n x (p * nbasis)} matrix of smoothed and centered
coefficients for the functional predictors.}
\item{eigvalues}{A vector of eigenvalues obtained from the eigenvalue
decomposition of the matrix of interest that is calculated during the
dimension reduction process.}
\item{phi}{A \code{(p * nbasis) x (p * nbasis)} matrix of eigenvectors
resulting from the eigenvalue decomposition of the matrix of interest
that is calculated during the dimension reduction process.}
\item{gx}{The \code{(p * nbasis) x (p * nbasis)} block diagonal Gram
matrix of the B-spline basis functions.}
\item{betas}{The \code{(p * nbasis) x (p * nbasis)} matrix of the
coordinates of \eqn{\beta_1, \dots, \beta_d} resulting from the
coordinate representation on the B-spline basis functions.}
\item{sufpred}{A \code{n x (p * nbasis)} matrix of estimated sufficient
predictors of the functional central subspace.}
}
\description{
\code{mfsir} performs dimension reduction for functional data and
provides the new estimated predictors.
}
\details{
This function implements functional sliced inverse regression (FSIR),
introduced by Ferré and Yao (2003), for scalar-on-function
regression problems.  The method seeks to identify a low-dimensional
subspace that captures the relationship between a scalar response and
functional predictors.  The authors proved that \eqn{E(X|Y) - E(X)}
belongs to \eqn{\Sigma_{XX} S_{Y|X}}, where \eqn{S_{Y|X}} denotes the
functional central subspace.

For \eqn{i=1, \dots, p}, let \eqn{\mathcal{H}_i} be a separable
Hilbert space of real-valued functions on \eqn{T}, a bounded closed interval
in \eqn{\mathbb{R}}.  Let \eqn{Y: \Omega \rightarrow \mathbb{R}} a univariate
response and \eqn{X=(X^1, \dots, X^p): \Omega \rightarrow
\bigoplus_{i=1}^p \mathcal{H}_i} a random element.  Dimension reduction
techniques aim at finding functions \eqn{\beta_1, \dots, \beta_d} in
\eqn{\bigoplus_{i=1}^p \mathcal{H}_i}, such that
\deqn{Y = g(\langle \beta_1, X \rangle_{\bigoplus \mathcal{H}}, \dots,
\langle \beta_d, X \rangle_{\bigoplus \mathcal{H}}, \epsilon),} where
\eqn{g} is an arbitrary unknown function on \eqn{\mathbb{R}^{d+1}}, and
\eqn{\epsilon} is independent of \eqn{X}.  This implies that \eqn{Y} and
\eqn{X} are independent given \eqn{\langle \beta_1,
X \rangle_{\bigoplus \mathcal{H}}, \dots,
\langle \beta_d, X \rangle_{\bigoplus \mathcal{H}}} and that the \eqn{p}-
dimensional predictor \eqn{X} can be replaced with the \eqn{d}-dimensional
predictor \eqn{\langle \beta_1, X \rangle_{\bigoplus \mathcal{H}}, \dots,
\langle \beta_d, X \rangle_{\bigoplus \mathcal{H}}}.

The functions \eqn{\beta_1, \dots, \beta_d} are called the \emph{functional
dimension reduction directions} and the subspace spanned by \eqn{\beta_1,
\dots, \beta_d} is called the \emph{functional dimension reduction
subspace}.  The smallest functional dimension reduction subspace is called
the \emph{functional central subspace} and is denoted by \eqn{S_{Y|X}}.
}
\examples{
# set the parameters
n <- 100
p <- 5
nt <- 100
nbasis <- 4
H <- 10
tt <- seq(0, 1, length.out = nt)
eta.mat <- mvtnorm::rmvnorm(n, mean = rep(0, p * nbasis))
eta <- array(eta.mat, dim = c(nbasis, n, p))
# Generate the functional data
data <- fundata(n, p, nbasis, tt, 'fourier', eta)
Xc <- data$Xc
# Generate the model
error <- rnorm(n)
y <- 3 * data$mfpca.scores[, 1] + error
# Run mfsir
result <- mfsir(Xc, y, H, nbasis)
result$sufpred
# Plot the first sufficient predictor against the true one
plot(result$sufpred[, 1], data$mfpca.scores[, 1], xlab = 'First
    Sufficient Predictor', ylab = 'First true Predictor', pch = 20)
# Calculate the correlation between the estimated and true predictors
mcorr(result$sufpred[, 1], data$mfpca.scores[, 1])

}
\references{
Ferré, L, and Yao, F. (2003) Function Sliced Inverse Regression
Analysis. \emph{Statistics}, 37(6), 475-488.
}

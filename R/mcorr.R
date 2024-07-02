#' Multiple correlation
#'
#' \code{mcorr} computes the multiple correlation between two vectors or
#' matrices \code{u} and \code{v}.
#'
#' This function computes the multiple correlation between two vectors or
#' matrices \code{u} and \code{v} of the same dimension.  Let \eqn{C_{uu}},
#' \eqn{C_{uv}}, and \eqn{C_{vu}} represent the sample covariance matrices.
#' The multiple correlation between \code{u} and \code{v} is given by:
#' \deqn{
#' mcorr(u, v) = tr(C_{vv}^{-1/2} C_{vu} C_{uu}^{-1} C_{uv} C_{vv}^{-1/2}).
#' }
#' This number ranges from 0 to d, the dimension of \code{u} and \code{v},
#' where a number closer to d indicates stronger correlation.  This measure
#' was used, for example, in Li and Song (2022) and Solea et al. (2024).
#'
#' @param u,v Numeric vectors or matrices.
#' @return \code{mcorr} computes the multiple correlation between \code{u}
#'    and \code{v}.
#'
#' @references Li, B, and Song, J. (2022). Dimension reduction for functional
#'    data based on weak conditional moments. \emph{The Annals of Statistics}
#'    50(1), 107--128.
#'
#' Solea, E., Christou, E., and Song, J. (2024). Robust Inverse Regression for
#'     Multivariate Elliptical Functional Data. \emph{Statistica Sinica}
#'     https://www3.stat.sinica.edu.tw/ss_newpaper/SS-2023-0341_na.pdf
#'
#' @include matpower.R
#' @examples
#' # Example 1
#' u <- matrix(rnorm(100), ncol = 2)
#' v <- u + 0.1 * matrix(rnorm(100), ncol = 2)
#' mcorr(u, v)
#'
#' # Example 2
#' u <- matrix(rnorm(100), ncol = 3)
#' v <- matrix(rnorm(100), ncol = 3)
#' mcorr(u, v)
#'
#' @export
mcorr <- function(u, v) {
  u <- as.matrix(u)
  v <- as.matrix(v)

  # Check if the number of rows for u and v agree
  if (nrow(u) != nrow(v)) {
    stop("u and v must have the same number of rows.")
  }

  # Check if u and v have the same dimension
  if (ncol(u) != ncol(v)) {
    stop("u and v must have the same number of columns.")
  }

  # Compute the correlation
  if(ncol(u) == 1) {
    result <- c(abs(stats::cor(u, v)))
    return(result)
  } else {
    suu.nhalf <- matpower(stats::var(u), -1 / 2)
    svv.inv <- solve(stats::var(v))
    result <- c(abs(sum(diag(suu.nhalf %*% stats::cov(u, v) %*% svv.inv %*%
                               stats::cov(v, u) %*% suu.nhalf))))
    return(result)
    }
}

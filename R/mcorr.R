#' Multiple correlation
#'
#' \code{mcorr} computes the multiple correlation between two random vectors
#' \code{u} and \code{v}.
#'
#' This function computes the multiple correlation between two random vectors
#' \code{u} and \code{v} of the same dimension. Let C_{uu}, C_{uv}, and C_{vu}
#' represent the sample covariance matrices, the multiple correlation between
#' \code{u} and \code{v} is
#' \deqn{
#' mcorr(u, v) = tr(C_{vv}^{-1/2} C_{vu} C_{uu}^{-1} C_{uv} C_{vv}^{-1/2}).
#' }
#' This was used, for example, in Li and Song (2022) and Solea et al. (2024).
#'
#' @param u A numeric vector or matrix. The first variable for the correlation
#'    measure.
#' @param v A numeric vector or matrix. The second variable for the correlation
#'    measure.
#' @return \code{mcorr} computes the multiple correlation between \code{u}
#'    and \code{v}.
#'
#' @references Li, B, and Song, J. (2022). Dimension reduction for functional
#'    data based on weak conditional moments. \emph{The Annals of Statistics}
#'    50(1), 107--128.
#'
#' @import stats
#' @include matpower.R
#' @examples
#' set.seed(1234)
#' u <- matrix(rnorm(100), ncol=2)
#' v <- matrix(rnorm(100), ncol=2)
#' mcorr(u, v)
#'
#' @export
mcorr <- function(u, v) {
  u <- as.matrix(u)
  v <- as.matrix(v)

  # compatibility checks
  # checks if the number of rows for u and v agree
  if (dim(u)[1] != dim(v)[1]) {
    stop("u and v must have the same number of rows.")
  }

  # checks if u and v have the same dimension
  if (dim(u)[2] != dim(v)[2]) {
    stop("u and v must have the same dimension.")
  }

  if(dim(u)[2] == 1) {
    result <- c(abs(cor(u, v)))
    return(result)
  }

  if(dim(u)[2] != 1) {
    suu.nhalf <- matpower(var(u), -1/2)
    svv.inv <- solve(var(v))
    result <- c(abs(sum(diag(suu.nhalf %*% cov(u, v) %*% svv.inv %*%
                               cov(v, u) %*% suu.nhalf))))
    return(result)
    }
}

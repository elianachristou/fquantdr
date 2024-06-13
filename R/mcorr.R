#' multiple correlation
#'
#' \code{mcorr} calculates the multiple correlation between two random vectors
#' \code{u} and \code{v}
#'
#'



mcorr <- function(u, v) {
  u <- as.matrix(u)
  v <- as.matrix(v)
  if(dim(u)[2] == 1) return(c(abs(cor(u, v))))
  if(dim(u)[2] != 1) {
    suu.nhalf <- matpower(var(u), -1/2)
    svv.inv <- solve(var(v))
    return(c(abs(sum(diag(suu.nhalf %*% cov(u, v) %*% svv.inv %*% cov(v, u) %*% suu.nhalf)))))}
}

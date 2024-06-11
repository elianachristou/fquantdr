#' Discretize Response
#'
#' \code{discretize} converts the response into a discrete form by creating slices and assigning each part of the response to the slice it is in.
#'
#' This function converts a vector response into a vector of the slices each value in the response is a part of. This is important for the Functional Sliced Inverse Regression method which relies on using these slices for reducing the dimension of data.
#'
#' @param y A vector representing the response variable
#' @param yunit A vector defining the slices used
#'
#' @return A vector defining the slices each item in the response corresponds to
#'
#' @noRd
#' @examples
#' y -> c(1, 2, 3, 2)
#' yunit -> unique(y)
#'
discretize <- function(y, yunit) {
  n <- length(y)
  #Add small amount of noise to y
  y <- y + .00001 * mean(y) * rnorm(n)
  nsli <- length(yunit)
  #order y values in ascending order
  yord <- y[order(y)]
  n <- length(y)
  nwith <- floor(n / nsli)
  #instantiate vector of division points between slices
  divpt <- rep(0, nsli - 1)
  #set each division point to the values of yord that
  #represent the boundaries between slices
  for(i in 1:(nsli - 1)) {
    divpt[i] <- yord[i * nwith + 1] }

  y1 <- rep(0, n)
  #Assign slice labels to the upper boundary slice
  y1[y >= divpt[nsli - 1]] <- nsli
  #Assign slice labels to the lower boundary slice
  y1[y < divpt[1]] <- 1
  #Assign slices labels to intermediate slices
  for(i in 2:(nsli - 1)) {
    y1[(y >= divpt[i - 1]) & (y < divpt[i])] <- i }
  #return discretized y
  return(y1)
}

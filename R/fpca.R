#' Functional Principal Component Analysis
#'
#' \code{fpca} performs functional principal component analysis (FPCA) on functional
#' data.
#'
#' This function performs FPCA on the given functional data and returns the
#' eigenvalues and eigenfunctions of the covariance matrix.
#'
#' @param ftn A list containing the coefficients of the functional object and
#'    the basis
#' @param basisname A string specifying the type of basis function to use for FPCA
#'
#' @return \code{fpca} returns:
#'    \itemize{
#'        \item \code{pred}: The principal component scores matrix
#'        \item \code{eval}: The eigenvalues
#'        \item \code{mat}: The matrix of eigenfunctions
#'    }

#' Functional data
#'
#' \code{fundata} generates functional data based on basis functions.
#'
#' This function generates functional data for a given number of observations
#' ('n'), functional predictors ('p'), basis functions ('nbasis'), and time
#' points ('time').  It utilizes a specified coefficient matrix ('eta') to
#' create the functional predictors.  Using these basis functions, it computes
#' the original functional predictors ('x') by multiplying the coefficients
#' with the basis functions. The function then centers these functional
#' predictors ('xc') by subtracting the mean across observations for each time
#' point.  The output is a list containing both the original and centered
#' functional predictors.
#'
#' @param n The number of observations, i.e., sample size.
#' @param p The number of functional predictors.
#' @param nbasis The number of Fourier basis functions.
#' @param time A vector of time points at which the functional data is
#'     evaluated.
#' @param eta A matrix of coefficients with dimensions \code{n} by
#'     \code{p * nbasis}.
#'
#' @return \code{fundata} generates functional data based on Fourier
#'    basis functions and returns:
#'    \itemize{
#'        \item \code{x}: The original functional predictors, a
#'            \code{n * nt * p} array, where \code{nt} denotes the number
#'            of time points.
#'        \item \code{xc}: The centered functional predictors, a
#'            \code{n * nt * p} array, where \code{nt} denotes the number
#'            of time points.
#'    }
#'

# ------------------------
# Correlator Documentation
# ------------------------
#' @title Correlation Structure
#'
#' @description Creates a correlation structure u(x) for an emulator.
#'
#'     Requires two functions - one that returns the correlation between two points,
#'     and one which calculates the expectation of a given point under this correlation function.
#'     Both can be constants, but must be expressed as a function. If a nugget is required, then
#'     the form of the correlator is changed: without a nugget the equation for Cov[u(x),u(x')]
#'     is simply \code{cov}; with a nugget this changes to \code{(1-delta)*cov + delta*cov(0,0)*I}
#'     with \code{I} being an indicator function.
#'
#' @name Correlator
#'
#' @section Constructor: \code{Correlator$new(cov = NULL, mu = NULL, delta = NULL)}
#'
#' @section Arguments:
#'  \code{cov} A covariance function, which must take two arguments.
#'
#'  \code{mu} An expectation function.
#'
#'  \code{delta} The size of the 'nugget' to be added, as a numeric in the range [0,1).
#'
#' @section Constructor Details: The constructor must have both a covariance function and
#' an expectation function of the correct type. Errors will be thrown if \code{cov} does not
#' take two arguments, and similarly if \code{mu} does not take one. The nugget is optional
#' (but may be useful!)
#'
#' @section Public Methods:
#'  \code{$get_exp(x)} Returns the expectation of the correlator at a point x.
#'
#'  \code{$get_cov(x, xp = NULL)} Returns the covariance between two points, or the variance
#'  at a point if only one argument is supplied.
#'
#' @export
#'
#' @examples
#'     expectation_function <- function(x) 1
#'     covariance_function <- function(x, xp) 0.36*exp_sq(x, xp, 0.1)
#'     u <- Correlator$new(covariance_function, expectation_function)
#'     u$get_exp(c(0,0.1,0.2))
#'     #> 1
#'     u$get_cov(c(0,0.1,0.3))
#'     #> 0.36
#'     u$get_cov(c(0,0.1,0.3),c(-0.1,0,0.2))
#'     #> 0.01792334
#'     # Now add a nugget term:
#'     delta = 0.08522
#'     unew <- Correlator$new(covariance_function, expectation_function, delta)
#'     unew$get_cov(c(0,0.1,0.3))
#'     #> 0.36 (as expected; nugget does not change variance)
#'     unew$get_cov(c(0,0.1,0.3), c(-0.1,0,0.2))
#'     #> 0.01639592
NULL
# ---------------------
# Correlator Definition
# ---------------------
Correlator <- R6::R6Class("Correlator", lock_objects = F)

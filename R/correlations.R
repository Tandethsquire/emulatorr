#' Exponential squared correlation function
#'
#' For points \code{x}, \code{xp} and a correlation length \code{theta}, gives the exponent
#' of the squared distance between \code{x} and \code{xp}, weighted by \code{theta} squared.
#'
#' @param x A numeric position vector
#' @param xp A numeric position vector
#' @param theta A numeric correlation length
#'
#' @return The exponental-squared correlation between x and xp.
#' @export
#' @examples
#' exp_sq(1,2,0.1)
#' #> 3.720076e-44
#' exp_sq(c(1,2,-1),c(1.5,2.9,-0.7),0.2)
#' #> 3.266131e-13
exp_sq <- function(x, xp, theta) {
  return(exp(-sum((x-xp)^2/theta^2)))
}

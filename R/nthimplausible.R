#' n-th Maximum Implausibility
#'
#' For a collection of emulators, it can be helpful to combine the implausibility
#' measures for a given observation. The maximum implausibility is, simply, the largest
#' implausibility value given by the emulators for each output; the 2nd maximum is the maximum
#' of the set without the maximum, and so on.
#'
#' @param emulators A set of \code{\link{Emulator}} objects.
#' @param x An input point
#' @param z The observed outputs, either as a numeric vector or as a collection of val, sigma pairs (see examples)
#' @param n The implausibility level to return. By default, the median implausibility is chosen
#' @param max_imp A maximum implausibility to consider: in most cases, it is useful to truncate the size of the I(x). Default: 20.
#'
#' @return The n-th maximum implausibility value.
#' @export
#'
#' @examples
#' ems <- emulator_from_data(GillespieSIR, output_names = c('nS', 'nI', 'nR'),
#'  ranges = list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05)),
#'  quadratic = TRUE)
#' targets <- list(
#'  list(val = 281, sigma = 10.43),
#'  list(val = 30, sigma = 11.16),
#'  list(val = 689, sigma = 14.32)
#' )
#' nth_implausible(ems, data.frame(aSI = 0.4, aIR = 0.25, aSR = 0.025), targets)
#' grid <- expand.grid(
#'  aSI = seq(0.1, 0.8, length.out = 4),
#'  aIR = seq(0, 0.5, length.out = 4),
#'  aSR = seq(0, 0.05, length.out = 4)
#' )
#' nth_implausible(ems, grid, targets, n = 2)
nth_implausible <- function(emulators, x, z, n = 1, max_imp = 20) {
  if (is.numeric(z))
    output <- purrr::map(z, ~list(val=.x, sigma=0))
  else
    output <- z
  implausible_list <- t(apply(matrix(unlist(purrr::map(seq_along(emulators), ~emulators[[.x]]$implausibility(x, z[[.x]]))), ncol = length(emulators)), 1, sort))
  return(purrr::map_dbl(implausible_list[,length(implausible_list[1,])-n+1], ~min(.x, max_imp)))
}

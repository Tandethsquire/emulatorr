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
#'     out_vars <- c("nS", "nI", "nR")
#'     c_lengths <- c(0.1, 0.085, 0.075)
#'     ranges <- list(aSI = c(0.1,0.8), aIR = c(0,0.5), aSR = c(0,0.05))
#'     base_emulators <- emulator_from_data(GillespieSIR, names(ranges),
#'      out_vars, ranges = ranges, c_lengths = c_lengths)
#'     trained_emulators <- purrr::map(seq_along(base_emulators),
#'      ~base_emulators[[.x]]$bayes_adjust(GillespieSIR[,names(ranges)],
#'       GillespieSIR[,out_vars[[.x]]]))
#'    target_vals <- c(281, 30, 689)
#'    target_sigmas <- c(37.26, 11.16, 31.72)
#'    z_specs <- purrr::map2(target_vals, target_sigmas, ~list(val=.x, sigma=.y))
#'    nth_implausible(emulators = trained_emulators, x = c(0.4, 0.25, 0.025),
#'     z = z_specs, n = 2)
nth_implausible <- function(emulators, x, z, n = floor(length(emulators)/2), max_imp = 20) {
  if (is.numeric(z))
    output <- purrr::map(z, ~list(val=.x, sigma=0))
  else
    output <- z
  implausible_list <- sort(purrr::map_dbl(seq_along(emulators), ~min(emulators[[.x]]$implausibility(x, output[[.x]]),max_imp)))
  return(implausible_list[length(implausible_list)-n+1])
}

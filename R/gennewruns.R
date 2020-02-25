#' Generate Simulator Runs
#'
#' Given a set of trained emulators, finds the next set of points that will be informative
#' for the next wave of emulators.
#'
#' If the method is LHS, this creates a new training set using LHS, and
#' then finds the trace of the variance matrix of the emulators across these points (this is
#' broadly equivalent to using V-optimality). We repeat this for \code{n_runs}, and select
#' the configuration that minimises the mean of the variances across the emulators.
#' If observations are given, then these are used to ensure that the new sample points are
#' non-implausible according to the current emulators.
#'
#' Alternatively, we may use slice sampling: given the emulators and a known
#' non-implausible point, generate more by sampling uniformly over the minimum enclosing
#' hyperrectangle, and ensuring that the resulting point is implausible.
#'
#' For any sampling strategy, the parameters \code{emulators} and \code{ranges} must be
#' specified. If the method is 'slice', then the parameters \code{x}
#' and \code{z} are necessary. All other parameters are optional.
#'
#' @importFrom stats setNames runif
#'
#' @param emulators A list of \code{\link{Emulator}} objects, trained on the design points
#' @param ranges The ranges of the input parameters
#' @param n_points Optional. Specifies how many additional points are required. Default: 10*(number of emulators)
#' @param n_runs Optional. Number of sampling runs to perform. Default = 100
#' @param z Optional. If given, checks implausibility of sample points to restrict to only non-implausible points.
#' @param cutoff Optional. If z is given, this is the implausibility cutoff for the filtering
#' @param x Required if method is 'slice'. A starting point for the sampling
#' @param method 'lhs' or 'slice'
#'
#' @return A \code{data.frame} containing the set of new points to simulate at.
#' @export
#'
#' @examples
#'    ranges <- list(aSI = c(0.1,0.8), aIR = c(0,0.5), aSR = c(0,0.05))
#'    out_names <- c('nS', 'nI', 'nR')
#'    emulators <- emulator_from_data(input_data = GillespieSIR,
#'     input_names = names(ranges), output_names = c('nS','nI','nR'),
#'     ranges = ranges)
#'    trained_emulators <- purrr::map(seq_along(emulators),
#'     ~emulators[[.x]]$bayes_adjust(GillespieSIR[,names(ranges)],
#'     GillespieSIR[,out_names[[.x]]]))
#'    generate_new_runs(trained_emulators,
#'      ranges = ranges, n_points = 5, n_runs = 5)
#'    x_initial <- unlist(data.frame(aSI = 0.4, aIR = 0.25, aSR = 0.025))
#'    z <- list(
#'       list(val = 457, sigma = 15),
#'       list(val = 81, sigma = 7),
#'       list(val = 462, sigma = 16))
#'    generate_new_runs(trained_emulators, ranges,
#'       x = x_initial, z = z, method = "slice", n_points = 5)
generate_new_runs <- function(emulators, ranges, n_points = 10*length(ranges), n_runs = 20, z, cutoff = 3, x = NULL, method = "lhs") {
  if (method == "lhs")
  {
    current_trace = NULL
    out_points = NULL
    new_points <- t(apply(2*(lhs::optimumLHS(n_points*5, length(ranges))-1/2), 1, scale_input, ranges, FALSE))
    if (!missing(z)) new_points <- new_points[apply(new_points, 1, function(x) nth_implausible(emulators, x, z, 2, 20))<=cutoff,]
    if (length(new_points[,1]) < n_points) {
      message(cat("Only", length(new_points[,1]), "points generated."))
      return(setNames(data.frame(new_points), names(ranges)))
    }
    for (i in 1:n_runs)
    {
      new_points <- new_points[sample(seq_along(new_points[,1]), n_points),]
      measure <- mean(purrr::map_dbl(seq_along(emulators), ~sum(apply(new_points, 1, emulators[[.x]]$get_cov))))
      if (is.null(current_trace) || measure < current_trace) {
        out_points <- new_points
        current_trace <- measure
      }
    }
    return(setNames(data.frame(out_points), names(ranges)))
  }
  if (method == "slice") {
    out_points <- x
    x0 <- c(out_points, use.names = F)
    x_new = x0
    for (i in 1:n_points) {
      for (j in 1:length(ranges)) {
        xl = ranges[[j]][1]
        xr = ranges[[j]][2]
        repeat {
          x_new[j] <- runif(1, min = xl, max = xr)
          if (nth_implausible(emulators, x_new, z, n=2)<=cutoff) break
          else ifelse(x_new[j]<x0[j], xl <- x_new[j], xr <- x_new[j])
        }
      }
      out_points <- cbind(out_points, x_new)
      x0 <- x_new
    }
    return(setNames(data.frame(t(out_points), row.names = NULL), names(ranges)))
  }
}

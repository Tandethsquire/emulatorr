#' Generate Simulator Runs
#'
#' Given a set of trained emulators, finds the next set of points that will be informative
#' for the next wave of emulators. It first augments the current training set using LHS, and
#' then finds the trace of the variance matrix of the emulators across these points (this is
#' broadly equivalent to using V-optimality). We repeat this for \code{n_runs}, and select
#' the configuration that minimises the mean of the variances across the emulators.
#'
#' @param old_points A \code{data.frame} containing the original design points.
#' @param emulators A list of \code{\link{Emulator}} objects, trained on the design points
#' @param ranges The ranges of the input parameters
#' @param n_points Optional. Specifies how many additional points are required. Default: 10*(number of emulators)
#' @param n_runs Optional. Number of sampling runs to perform. Default = 100
#'
#' @return A \code{data.frame} containing the set of new points to simulate at.
#' @export
#'
#' @examples
#'    emulators <- emulator_from_data(input_data = GillespieSIR,
#'     input_names = c('aSI','aIR','aSR'), output_names = c('nS','nI','nR'),
#'     c_lengths = c(0.1, 0.085, 0.075))
#'    trained_emulators <- purrr::map(seq_along(emulators),
#'     ~emulators[[.x]]$bayes_adjust(GillespieSIR[,c('aSI','aIR','aSR')],
#'     GillespieSIR[,c('nS','nI','nR')[[.x]]]))
#'    gen_new_runs(GillespieSIR[,c('aSI','aIR','aSR')], trained_emulators,
#'     list(c(0.1, 0.8), c(0, 0.5), c(0, 0.05)), n_points = 5, n_runs = 10)
gen_new_runs <- function(old_points, emulators, ranges, n_points = 10*length(emulators), n_runs = 100) {
  range_lengths <- purrr::map_dbl(ranges, ~.x[2]-.x[1])
  orig_points <- as.matrix(old_points)
  current_trace = NULL
  out_points = NULL
  for (i in 1:n_runs) {
    combined_points <- lhs::augmentLHS(orig_points, n_points)
    new_points <- combined_points[(length(orig_points[,1])+1):length(combined_points[,1]),]
    new_points <- t(apply(new_points, 1, function(x) x*range_lengths + purrr::map_dbl(ranges, ~.x[1])))
    measure <- mean(purrr::map_dbl(seq_along(emulators), ~sum(apply(new_points, 1, emulators[[.x]]$get_var))))
    if (is.null(current_trace) || measure < current_trace) {
      out_points <- new_points
      current_trace <- measure
    }
  }
  return(data.frame(out_points))
}

#' Emulator Standard Errors
#'
#' Finds and plots emulator standard errors.
#' For an emulator of a simulator function \code{f(x)}, and a validation data set \code{X},
#' finds the standard errors in the form \code{(f(x)-E[f(x)])/sqrt(Var[f(x)])},
#' where \code{E[f(x)]} is the emulator expectation, and \code{Var[f(x)]} is the emulator variance,
#' at each point \code{x} in \code{X}.
#'
#' @importFrom graphics hist
#'
#' @param emulator An \code{\link{Emulator}} object.
#' @param input_points A set of validation points.
#' @param output_points The outputs, \code{f(x)}, from the simulator.
#' @param output_name Optional. A name for the output.
#'
#' @return A list of standard errors.
#' @export
#'
#' @examples
#'     in_vars <- c("aSI", "aIR", "aSR")
#'     out_vars <- c("nS", "nI", "nR")
#'     c_lengths <- c(0.1, 0.085, 0.075)
#'     base_emulators <- emulator_from_data(GillespieSIR, in_vars,
#'      out_vars, c_lengths = c_lengths)
#'     trained_emulators <- purrr::map(seq_along(base_emulators),
#'      ~base_emulators[[.x]]$bayes_adjust(GillespieSIR[,in_vars], GillespieSIR[,out_vars[[.x]]]))
#'     standard_errors(trained_emulators[[1]], GillespieValidation[,in_vars],
#'      GillespieValidation[,'nS'], "nS")
standard_errors <- function(emulator, input_points, output_points, output_name) {
  errors <- (apply(input_points, 1, function(x) emulator$get_exp(x))-output_points)/apply(input_points, 1, function(x) sqrt(emulator$get_var(x)))
  if (missing(output_name))
    hist(errors, main = "Standard Errors")
  else
    hist(errors, xlab="Standard Error", main = paste("Standard Errors for Output:",output_name))
  return(errors)
}

#' Emulator Diagnostic Plot
#'
#' Produces a diagnostic plot of emulator output.
#' The emulator output \code{E[f(x)]} is plotted against the simulator output \code{f(x)},
#' with error bars given by the emulator standard deviation \code{sqrt(3*Var[f(x)])},
#' for each point \code{x} in a validation set \code{X}. Points whose emulator expectation
#' lies outside 3-sigma of the simulator output are shown in red, and those input points
#' are returned.
#'
#' @importFrom graphics abline arrows plot
#'
#' @param emulator An \code{\link{Emulator}} object.
#' @param input_points A set of validation points.
#' @param output_points The outputs, \code{f(x)}, from the simulator.
#' @param range The expected range of the output.
#' @param sd Numeric: the allowed number of standard deviations (default: 3).
#' @param output_name Optional. A name for the output.
#'
#' @return A list of points whose emulator value is outside the allowed standard deviation.
#' @export
#'
#' @examples
#'     in_vars <- c("aSI", "aIR", "aSR")
#'     out_vars <- c("nS", "nI", "nR")
#'     c_lengths <- c(0.1, 0.085, 0.075)
#'     base_emulators <- emulator_from_data(GillespieSIR, in_vars,
#'      out_vars, c_lengths = c_lengths)
#'     trained_emulators <- purrr::map(seq_along(base_emulators),
#'      ~base_emulators[[.x]]$bayes_adjust(GillespieSIR[,in_vars], GillespieSIR[,out_vars[[.x]]]))
#'     comparison_diagnostics(trained_emulators[[1]], GillespieSIR[,in_vars],
#'      GillespieSIR[,'nS'], list(c(0,1000),c(0,1000)))
comparison_diagnostics <- function(emulator, input_points, output_points, range, sd=3, output_name) {
  emulator_outputs <- apply(input_points, 1, emulator$get_exp)
  emulator_uncertainty <- apply(input_points, 1, function(x) sd*sqrt(emulator$get_var(x)))
  emulator_invalid <- purrr::map2_lgl(output_points>emulator_outputs+emulator_uncertainty,
                                      output_points<emulator_outputs-emulator_uncertainty,
                                      ~.x|.y)
  if (missing(output_name))
    title_str <- "Simulator/Emulator Comparison"
  else
    title_str <- paste("Simulator/Emulator Comparison for Output:", output_name)
  plot(output_points, emulator_outputs, pch=16, col=ifelse(emulator_invalid, 'red', 'black'),
         xlim = range[[1]], ylim=range[[2]], xlab='f(x)', ylab='E[f(x)]',
         panel.first = c(abline(a=0, b=1, col = 'green')),
         main = title_str)
  for (i in seq_along(input_points[,1]))
    arrows(x0 = output_points[[i]], y0 = emulator_outputs[[i]] - emulator_uncertainty[[i]],
           x1 = output_points[[i]], y1 = emulator_outputs[[i]] + emulator_uncertainty[[i]],
           col = ifelse(emulator_invalid[[i]], 'red', 'blue'),
           code = 3, angle = 90, length = .1)
  which_invalid <- input_points[emulator_invalid,]
  return(which_invalid)
}

#' Classification diagnostics
#'
#' Checks for emulator misclassifications.
#' Both the emulator implausibility and the simulator implausibility are computed, and
#' plotted against one another. Points for which the emulator implausibility is outside the desired
#' cut-off but for which the simulator implausibility is not are misclassification points,
#' and are highlighted in red.
#'
#' @importFrom graphics abline plot
#'
#' @param emulator An \code{\link{Emulator}} object.
#' @param input_points A set of validation points.
#' @param output_points The outputs, \code{f(x)}, from the simulator.
#' @param z The observation to test implausibility against. Either as a single \code{numeric}, or as \code{list(val=numeric, sigma=numeric)}.
#' @param output_name Optional. A name for the output.
#' @param cutoff Optional. The cut-off for the implausibility measure.
#'
#' @return The set of points misclassified by the emulator.
#' @export
#'
#' @examples
#'     in_vars <- c("aSI", "aIR", "aSR")
#'     out_vars <- c("nS", "nI", "nR")
#'     c_lengths <- c(0.1, 0.085, 0.075)
#'     base_emulators <- emulator_from_data(GillespieSIR, in_vars,
#'      out_vars, c_lengths = c_lengths)
#'     trained_emulators <- purrr::map(seq_along(base_emulators),
#'      ~base_emulators[[.x]]$bayes_adjust(GillespieSIR[,in_vars], GillespieSIR[,out_vars[[.x]]]))
#'     target_value <- list(val = 281, sigma = 37.26)
#'     classification_error(emulator = trained_emulators[[1]],
#'     input_points = GillespieValidation[,in_vars], output_points <- GillespieValidation[,'nS'],
#'     z = target_value, output_name = 'nS')
classification_error <- function(emulator, input_points, output_points, z, output_name, cutoff=3)
{
  if (is.numeric(z))
    output <- list(val = z, sigma = 0)
  else
    output <- z
  emulator_implausibility <- apply(input_points, 1, function(x) emulator$implausibility(x, output))
  simulator_implausibility <- purrr::map_dbl(output_points, ~sqrt((output$val-.x)^2/output$sigma^2))
  misclass <- purrr::map2_lgl(emulator_implausibility > cutoff, simulator_implausibility <= cutoff, ~.x&.y)
  if (missing(output_name))
    title_str <- "Emulator/Simulator Classification"
  else
    title_str <- paste("Emulator/Simulator Classification for Output:", output_name)
  plot(emulator_implausibility, simulator_implausibility, pch = 16,
       col = ifelse(misclass, 'red', 'black'), xlab = "Emulator Implausibility", ylab = "Simulator Implausibility",
       main = title_str, panel.first = c(abline(h=cutoff), abline(v=cutoff)))
  which_misclass <- input_points[misclass,]
  return(which_misclass)
}

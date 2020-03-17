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
#' @param plt Should a plot be shown (default: T).
#'
#' @return A list of standard errors.
#' @export
#'
#' @examples
#' em <- emulator_from_data(GillespieSIR, output_names = c('nS', 'nI', 'nR'),
#'  ranges = list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05)),
#'  quadratic = TRUE)[[1]]
#' standard_errors(em, GillespieValidation[,1:3], GillespieValidation[,'nS'], 'nS')
#' #> (0.7864384, 0.01426296, 0.001072935)
standard_errors <- function(emulator, input_points, output_points, output_name, plt = T) {
  errors <- (emulator$get_exp(input_points)-output_points)/sqrt(emulator$get_cov(input_points))
  if (plt) {
    if (missing(output_name))
      hist(errors, main = "Standard Errors")
    else
      hist(errors, xlab="Standard Error", main = paste("Standard Errors for Output:",output_name))
  }
  return(input_points[abs(errors)>3,])
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
#' @param sd Numeric: the allowed number of standard deviations (default: 3).
#' @param output_name Optional. A name for the output.
#' @param plt Should a plot be shown (default: T).
#'
#' @return A list of points whose emulator value is outside the allowed standard deviation.
#' @export
#'
#' @examples
#' em <- emulator_from_data(GillespieSIR, output_names = c('nS', 'nI', 'nR'),
#'  ranges = list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05)),
#'  quadratic = TRUE)[[1]]
#' comparison_diagnostics(em, GillespieValidation[,1:3], GillespieValidation[,'nS'])
#' #> (0.7864384, 0.01426296, 0.001072935)
comparison_diagnostics <- function(emulator, input_points, output_points, sd=3, output_name, plt=T) {
  emulator_outputs <- emulator$get_exp(input_points)
  emulator_uncertainty <- sd*sqrt(emulator$get_cov(input_points))
  em_ranges <- range(c(emulator_outputs + emulator_uncertainty, emulator_outputs - emulator_uncertainty))
  emulator_invalid <- purrr::map2_lgl(output_points>emulator_outputs+emulator_uncertainty,
                                      output_points<emulator_outputs-emulator_uncertainty,
                                      ~.x|.y)
  if (missing(output_name))
    title_str <- "Simulator/Emulator Comparison"
  else
    title_str <- paste("Simulator/Emulator Comparison for Output:", output_name)
  if (plt) {
    plot(output_points, emulator_outputs, pch=16, col=ifelse(emulator_invalid, 'red', 'black'),
           xlim = range(output_points), ylim=range(em_ranges), xlab='f(x)', ylab='E[f(x)]',
           panel.first = c(abline(a=0, b=1, col = 'green')),
           main = title_str)
    for (i in seq_along(input_points[,1]))
      arrows(x0 = output_points[[i]], y0 = emulator_outputs[[i]] - emulator_uncertainty[[i]],
             x1 = output_points[[i]], y1 = emulator_outputs[[i]] + emulator_uncertainty[[i]],
             col = ifelse(emulator_invalid[[i]], 'red', 'blue'),
             code = 3, angle = 90, length = .1)
  }
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
#' @param plt Should a plot be shown (default: T).
#'
#' @return The set of points misclassified by the emulator.
#' @export
#'
#' @examples
#' em <- emulator_from_data(GillespieSIR, output_names = c('nS', 'nI', 'nR'),
#'  ranges = list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05)),
#'  quadratic = TRUE)[[1]]
#' target_value <- list(val = 281, sigma = 37.26)
#' classification_error(em, GillespieValidation[,1:3], GillespieValidation[,'nS'], target_value)
#' #> data.frame containing 0 points
classification_error <- function(emulator, input_points, output_points, z, output_name, cutoff=3, plt=T)
{
  if (is.numeric(z))
    output <- list(val = z, sigma = 0.001)
  else
    output <- z
  emulator_implausibility <- emulator$implausibility(input_points, output)
  simulator_implausibility <- purrr::map_dbl(output_points, ~sqrt((output$val-.x)^2/output$sigma^2))
  misclass <- purrr::map2_lgl(emulator_implausibility > cutoff, simulator_implausibility <= cutoff, ~.x&.y)
  if (missing(output_name))
    title_str <- "Emulator/Simulator Classification"
  else
    title_str <- paste("Emulator/Simulator Classification for Output:", output_name)
  if (plt) {
    plot(emulator_implausibility, simulator_implausibility, pch = 16,
         col = ifelse(misclass, 'red', 'black'), xlab = "Emulator Implausibility", ylab = "Simulator Implausibility",
         main = title_str, panel.first = c(abline(h=cutoff), abline(v=cutoff)))
  }
  which_misclass <- input_points[misclass,]
  return(which_misclass)
}

#' Validation Set Plotting
#'
#' Plots each of the emulator outputs against each input. For each input parameter, the
#' emulator expectation is plotted for each output. These plots are presented as a set, to
#' better identify trends and dependencies in the emulators outputs.
#'
#' @importFrom graphics plot par
#' @importFrom stats setNames
#'
#' @param emulators A set of \code{\link{Emulator}} objects.
#' @param input_points The validation input points on which to evaluate the emulators
#' @param output_names The names of the output parameters.
#'
#' @return NULL
#' @export
#' @examples
#' ems <- emulator_from_data(GillespieSIR, output_names = c('nS', 'nI', 'nR'),
#'  ranges = list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05)),
#'  quadratic = TRUE)
#' validation_plots(ems, GillespieValidation[,1:3], c('nS', 'nI', 'nR'))
validation_plots <- function(emulators, input_points, output_names) {
  input_names <- names(emulators[[1]]$ranges)
  out_data <- setNames(as.data.frame(cbind(input_points, purrr::map(emulators, ~.x$get_exp(input_points)))), c(input_names, output_names))
  op <- par(mfrow = c(length(output_names), length(input_names)))
  for (i in 1:length(emulators)) {
    for (j in 1:length(input_names)) {
      plot(out_data[,input_names[[j]]], out_data[,output_names[[i]]], pch=16, xlab = input_names[[j]], ylab=output_names[[i]])
    }
  }
  par(op)
}

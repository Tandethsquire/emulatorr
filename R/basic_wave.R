#' History Match
#'
#' Performs a full wave of emulation and history matching from data.
#' Given simulator runs (split into training and validation data), the target values,
#' and the identification of outputs to emulate, the function generates trained emulators,
#' tests them with emulator diagnostics (removing any emulators whose outputs cannot be
#' well emulated from the data), and finally generates a new sample of points to be entered
#' into the simulator.
#'
#' Necessary parameters to be passed are the input data, the validation data, the ranges of
#' inputs, and the observation values for each output. If any specifications are to be passed
#' directly to the emulator construction, then they should be given as additional parameters
#' (see \code{\link{emulator_from_data}} to see the options).
#'
#' If the wave to be emulated is the first wave, then a set of preliminary ('wave 0') emulators
#' are fitted to the data. The proper set of emulators are then trained using Bayes Linear
#' adjustment. On subsequent waves, the preliminary emulators should be passed to the function
#' as the argument \code{previous_wave}.
#'
#' The output consists of a list of four items: the preliminary emulators \code{base_emulators},
#' the trained emulators \code{emulators}, the next points to be put into the simulator
#' \code{next_sample}, and the minimum enclosing hyperrectangle for the non-implausible
#' region, given as ranges \code{new_ranges}.
#'
#' @param input_data The set of training points
#' @param validation_data The set of points to use in validation
#' @param ranges The ranges of the inputs, as a named list.
#' @param output_names The names of the outputs to emulate.
#' @param targets The observations, given in the usual \code{(val, sigma)} form
#' @param n_points The number of points to evaluate the parameters on.
#' @param previous_wave The preliminary emulators for the set of waves, if they exist
#' @param sample_method The method to be used to find new points (see \code{\link{generate_new_runs}})
#' @param ... Any optional parameters to pass to \code{\link{emulator_from_data}}
#'
#' @return A list of base emulators, trained emulators for this wave, new sample points, and new ranges.
#' @export
#'
#' @examples
#'  #ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#'  #outputs <- c('nS','nI','nR')
#'  #targets <- list(
#'  # list(val = 281, sigma = 10.43),
#'  # list(val = 30, sigma = 11.16),
#'  # list(val = 689, sigma = 14.32)
#'  #)
#'  #wave1 <- full_wave(GillespieSIR, GillespieValidation, ranges, outputs, targets,
#'  # n_points = 30, deltas = rep(0.1, 3), quadratic = TRUE)

full_wave <- function(input_data, validation_data, ranges, output_names, targets, n_points = 40, previous_wave = NULL, sample_method = 'lhs', ...) {
  targets <- setNames(targets, output_names)
  if (is.null(previous_wave))
    base_emulators <- emulator_from_data(input_data, output_names, ranges, ...)
  else
    base_emulators <- previous_wave$base_emulators
  trained_emulators <- setNames(purrr::map(seq_along(base_emulators), ~base_emulators[[.x]]$adjust(input_data, output_names[[.x]])), output_names)
  cat("Running diagnostics...\n")
  diaglist <- list()
  for (i in 1:length(trained_emulators)) {
    diagdata <- comparison_diagnostics(trained_emulators[[i]], validation_data[,names(ranges)], validation_data[,output_names[[i]]], plt = F)
    diagdata <- rbind(diagdata, classification_error(trained_emulators[[i]], validation_data[,names(ranges)], validation_data[,output_names[[i]]], z = targets[[i]], plt = F))
    diagdata <- rbind(diagdata, standard_errors(trained_emulators[[i]], validation_data[,names(ranges)], validation_data[,output_names[[i]]], plt = F))
    diaglist[[output_names[i]]] = unique(diagdata)
  }
  for (i in 1:length(output_names)) {
    if (length(diaglist[[output_names[[i]]]][,1]) > 0.2*length(input_data[,1])) {
      cat("Cannot reliably emulate output ",output_names[[i]],".\n", sep="")
      trained_emulators[[output_names[[i]]]] <- NULL
      targets[[output_names[[i]]]] <- NULL
    }
  }
  if (length(trained_emulators) < length(output_names)/2) stop("Not enough outputs can be emulated.")
  cat("Completed diagnostics. Finding non-implausible region...\n")
  makeGrid <- function(ranges, npoints) {
    seqs <- purrr::map(ranges, ~seq(.x[[1]], .x[[2]], length.out = npoints))
    return(setNames(expand.grid(seqs), names(ranges)))
  }
  eval_grid <- makeGrid(ranges, n_points)
  imps <- nth_implausible(trained_emulators, eval_grid, targets)
  imp_data <- setNames(data.frame(cbind(eval_grid, imps)), c(names(ranges), "I"))
  p_set <- imp_data[imp_data$I<=3,]
  new_ranges <- lapply(p_set[,names(ranges)], function(x) c(min(x), max(x)))
  cat("Generating new sample points...\n")
  if (sample_method == 'lhs')
    new_points <- generate_new_runs(trained_emulators, new_ranges, z = targets)
  else if (sample_method == 'slice') {
    sample_point <- unlist(p_set[sample(seq_along(p_set), 1), names(ranges)], use.names = F)
    new_points <- generate_new_runs(trained_emulators, new_ranges, z = targets, method = 'slice', x = sample_point)
  }
  else if (sample_method == 'optical')
    new_points <- generate_new_runs(trained_emulators, ranges, z = targets, method = 'optical', plausible_set = p_set)
  else
    stop("Sampling method not recognised.")
  return(list(base_emulators = base_emulators, emulators = trained_emulators, next_sample = new_points, new_ranges = new_ranges))
}

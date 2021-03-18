#' Emulator Diagnostics
#'
#' Plots standard diagnostics for emulators.
#'
#' These diagnostics are based on having two datasets: a training set and a
#' validation set. The emulators will have been trained on the training set,
#' and the validation set is passed to the functions in this wrapper.
#'
#' The current options for diagnostics (with the codes for which_diag) are:
#'
#'   Standard Errors (se)
#'
#'   Comparison Diagnostics (cd)
#'
#'   Classification Error (ce)
#'
#'   All of the above (all)
#'
#' For details on each of these, see the help files for \code{\link{standard_errors}},
#' \code{\link{comparison_diagnostics}} and \code{\link{classification_error}} respectively.
#'
#' @importFrom purrr %>%
#'
#' @param emulators A list of emulators on which to perform diagnostics
#' @param validation_points The validation set
#' @param output_names The list of outputs to perform diagnostics on
#' @param targets If required, the list of observations for the outputs
#' @param which_diag Which diagnostics should be performed?
#' @param ... Any additional parameters to pass to the diagnostic tests.
#'
#' @return A data.frame containing the points that failed one or more diagnostic tests.
#'
#' @export
#'
#' @examples
#' output_names <- c('nS','nI','nR')
#' ems <- emulator_from_data(GillespieSIR, output_names,
#'  ranges = list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05)),
#'  quadratic = TRUE)
#' targets <- list(
#'  list(val = 281, sigma = 10.43),
#'  list(val = 30, sigma = 11.16),
#'  list(val = 689, sigma = 14.32)
#' )
#' validation_diagnostics(ems, GillespieValidation, output_names, targets = targets)
#' validation_diagnostics(ems, GillespieValidation, output_names, c('se','cd'))
#' validation_diagnostics(ems[1:2], GillespieValidation, output_names[1:2], 'ce', targets[1:2])
#' validation_diagnostics(ems, GillespieValidation, output_names,
#'  targets = targets, sd = 1, cutoff = 4)
#'
validation_diagnostics <- function(emulators, validation_points, output_names, which_diag = 'all', targets = NULL, ...) {
  if (!is.null(targets) && is.null(names(targets))) names(targets) <- output_names
  in_points <- validation_points[,names(emulators[[1]]$ranges)]
  out_points <- validation_points[,output_names]
  fail_point_list <- list()
  cl <- TRUE
  if (length(which_diag) == 1 && which_diag == 'all') actual_diag <- c('se','cd','ce')
  else {
    actual_diag <- which_diag[which_diag %in% c('se', 'cd', 'ce')]
    if (length(which_diag) != length(actual_diag)) warning(paste("Unrecognised diagnostics:",paste0(which_diag[!which_diag %in% c('se','cd','ce')], collapse = ", "),"\n\tValid diagnostic labels are cd, ce, se or all."))
  }
  mf <- length(actual_diag)
  if (purrr::some(actual_diag, ~. == 'ce') && is.null(targets)) {
    warning("No observations provided; cannot perform classification diagnostics.")
    cl <- FALSE
    mf <- mf - 1
  }
  op <- par(mfrow = c(3, mf))
  for (i in 1:length(emulators))
  {
    if (purrr::some(actual_diag, ~. == 'cd')) fail_point_list[[length(fail_point_list)+1]] <- comparison_diagnostics(emulators[[i]], in_points, out_points[[i]], output_names[[i]], targets = targets, ...)
    if (purrr::some(actual_diag, ~. == 'se')) fail_point_list[[length(fail_point_list)+1]] <- standard_errors(emulators[[i]], in_points, out_points[[i]], output_names[[i]], targets = targets, ...)
    if (purrr::some(actual_diag, ~. == 'ce') && cl) fail_point_list[[length(fail_point_list)+1]] <- classification_error(emulators[[i]], in_points, out_points[[i]], targets[[i]], output_names[[i]], ...)
  }
  par(op)
  failed_points <- unique(data.frame(do.call('rbind', fail_point_list)) %>% setNames(names(in_points)))
  return(failed_points)
}

#' Emulator Standard Errors
#'
#' Finds and plots emulator standard errors.
#'
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
#' @param targets The output targets (to check if failing points are relevant). Default: NULL
#' @param ... Dummy parameters (for compatibility with diagnostic wrapper)
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
standard_errors <- function(emulator, input_points, output_points, output_name, plt = T, targets = NULL, ...) {
  errors <- (emulator$get_exp(input_points)-output_points)/sqrt(emulator$get_cov(input_points))
  if (plt) {
    if (missing(output_name))
      hist(errors, main = "Standard Errors")
    else
      hist(errors, xlab="Standard Error", main = paste("Standard Errors for Output:",output_name))
  }
  emulator_invalid <- abs(errors) > 3
  if (!is.null(targets)) {
    this_target <- targets[[output_name]]
    point_invalid <- (output_points < this_target$val - 6*this_target$sigma) | (output_points > this_target$val + 6*this_target$sigma)
    emulator_invalid <- (!point_invalid) & emulator_invalid
  }
  return(input_points[emulator_invalid,])
}

#' Emulator Diagnostic Plot
#'
#' Produces a diagnostic plot of emulator output.
#'
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
#' @param targets The output targets (to check if failing points are relevant). Default: NULL
#' @param ... Dummy parameters (for compatibility with diagnostic wrapper)
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
comparison_diagnostics <- function(emulator, input_points, output_points, output_name, sd=3, plt=T, targets = NULL, ...) {
  emulator_outputs <- emulator$get_exp(input_points)
  emulator_uncertainty <- sd*sqrt(emulator$get_cov(input_points))
  em_ranges <- range(c(emulator_outputs + emulator_uncertainty, emulator_outputs - emulator_uncertainty))
  emulator_invalid <- purrr::map2_lgl(output_points>emulator_outputs+emulator_uncertainty,
                                      output_points<emulator_outputs-emulator_uncertainty,
                                      ~.x|.y)
  if (!is.null(targets)) {
    this_target <- targets[[output_name]]
    point_invalid <- (output_points < this_target$val - 2*sd*this_target$sigma) | (output_points > this_target$val + 2*sd*this_target$sigma)
    emulator_invalid <- (!point_invalid) & emulator_invalid
  }
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
#'
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
#' @param ... Dummy parameters (for compatibility with diagnostic wrapper)
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
classification_error <- function(emulator, input_points, output_points, z, output_name, cutoff=3, plt=T, ...)
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

#' Output Plotting
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
#' behaviour_plots(ems, GillespieValidation[,1:3], c('nS', 'nI', 'nR'))
behaviour_plots <- function(emulators, input_points, output_names) {
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

#' Expectation and Variance Visualisation
#'
#' Plots the points of a wave as a pairs plot, coloured by emulator expectation (lower) and
#' variance (upper) for each output.
#'
#' @import ggplot2
#' @importFrom GGally ggpairs wrap getPlot putPlot ggmatrix_gtable
#' @importFrom viridis scale_colour_viridis
#' @importFrom cowplot plot_grid get_legend
#'
#' @param ems The list of emulators
#' @param input_points The points on which to evaluate the emulators
#' @param output_names The list of outputs
#'
#' @return NULL
#'
#' @examples
#' ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#' targets <- list(
#'  list(val = 281, sigma = 10.43),
#'  list(val = 30, sigma = 11.16),
#'  list(val = 689, sigma = 14.32)
#' )
#' outputs <- c('nS','nI','nR')
#' ems <- emulator_from_data(GillespieSIR, outputs, ranges, deltas = rep(0.1, 3))
#' t_ems <- purrr::map(seq_along(ems), ~ems[[.]]$adjust(GillespieSIR, outputs[[.]]))
#' visualisation_plot(t_ems, GillespieSIR, outputs)
#'
#' @export
visualisation_plot <- function(ems, input_points, output_names) {
  inputs <- input_points[, names(ems[[1]]$ranges)]
  exp_vals <- setNames(data.frame(do.call('cbind', purrr::map(ems, ~.$get_exp(input_points)))), output_names)
  cov_vals <- setNames(data.frame(do.call('cbind', purrr::map(ems, ~.$get_cov(input_points)))), output_names)
  limfun <- function(data, mapping) {
    ggplot(data = data, mapping = mapping) +
      geom_point(cex = 2)
  }
  for (i in 1:length(output_names)) {
    temp_data <- setNames(data.frame(cbind(inputs, exp_vals[,output_names[i]], cov_vals[,output_names[i]])), c(names(inputs), 'exp', 'var'))
    g <- ggpairs(temp_data, columns = 1:length(inputs),
                 title = paste("Emulator Expectation and Variance:", output_names[i], sep = " "),
                 lower = list(continuous = wrap(limfun), mapping = aes(colour = temp_data[,'exp'])),
                 upper = list(continuous = wrap(limfun), mapping = aes(colour = temp_data[,'var'])),
                 diag = 'blank', progress = FALSE) +
      theme_minimal()
    for (j in 1:length(inputs)) {
      for (k in 1:length(inputs)) {
        plt <- getPlot(g, j, k)
        new_plt <- if(j < k) plt + scale_colour_viridis(option = 'plasma', name = "Var[f(x)]") else plt + scale_colour_viridis(option = 'magma', name = "E[f(x)]")
        g <- putPlot(g, new_plt, j, k)
      }
    }
    legends <- list(get_legend(g$plots[[2]]), rep(NULL, length(inputs)-2), get_legend(g$plots[[length(inputs)+1]]))
    print(suppressMessages(plot_grid(plotlist = list(ggmatrix_gtable(g), plot_grid(plotlist = legends, ncol = 1)), rel_widths = c(1, 0.1))))
  }
}

#' Space Removal
#'
#' Finds the proportion of space removed as a function of implausibility cut-off, and of structural
#' discrepancy, or changed variance.
#'
#' The reduction in space is found by evaluating over a p^d regular grid, where p is chosen by
#' \code{n_points} and d is the dimension of the input space. Larger values of \code{n_points}
#' will give a more accurate reflection of removed space, at high computational cost. For the
#' purpose of quick diagnostics, \code{n_points = 5} is acceptable.
#'
#' The parameter \code{modified} can take three options: 'disc' (default) corresponding to
#' model discrepancy, 'var' corresponding to emulator variance, or 'corr' corresponding to
#' correlation length. In the first case, the implausibilities are recalculated with the
#' original emulators; in the latter two cases, the emulators are re-trained with the new
#' specifications. For this reason, one should expect the 'var' and 'corr' options to be
#' more computationally intensive.
#'
#' The returned output is a \code{data.frame} consisting of the percentage of space
#' removed at each cutoff value, for each modifed value of the varied parameter. The main
#' result, however, is the accompanying plot of this information.
#'
#' @import ggplot2
#' @importFrom stats setNames
#' @importFrom viridis scale_colour_viridis
#' @importFrom plyr mutate
#'
#' @param emulators A set of \code{\link{Emulator}} objects.
#' @param validation_points The validation set used in this wave.
#' @param z The observations with which to match, as \code{list(val, sigma)} pairs.
#' @param n_points The number of points in each dimension of the grid.
#' @param u_mod The percentage differences in structural discrepancy to examine.
#' @param intervals The set of implausibility cut-offs to consider.
#' @param modified What parameter should be varied in the analysis?
#'
#' @return A ggplot object corresponding to the plot.
#' @export
#'
#' @examples
#' ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#' targets <- list(
#'  list(val = 281, sigma = 10.43),
#'  list(val = 30, sigma = 11.16),
#'  list(val = 689, sigma = 14.32)
#' )
#' outputs <- c('nS','nI','nR')
#' ems <- emulator_from_data(GillespieSIR, outputs, ranges, deltas = rep(0.1, 3), quadratic = TRUE)
#' t_ems <- purrr::map(seq_along(ems), ~ems[[.]]$adjust(GillespieSIR, outputs[[.]]))
#' names(t_ems) <- outputs
#' removal <- space_removed(ems, GillespieValidation, targets,
#'  n_points = 5, u_mod = seq(0.75, 1.25, by = 0.25), intervals = seq(2, 6, by = 0.1))
#'
space_removed <- function(emulators, validation_points, z, n_points = 10, u_mod = seq(0.8, 1.2, by = 0.1), intervals = seq(0, 10, length.out = 200), modified = 'disc') {
  value <- variable <- NULL
  in_names <- names(emulators[[1]]$ranges)
  z_vals <- purrr::map_dbl(z, ~.$val)
  z_sigs <- purrr::map_dbl(z, ~.$sigma)
  on_grid <- setNames(expand.grid(purrr::map(emulators[[1]]$ranges, ~seq(.[[1]], .[[2]], length.out = n_points))), in_names)
  imp_array <- array(0, dim = c(length(intervals), length(u_mod)))
  if (!(modified %in% c('disc', 'var', 'corr'))) {
    warning("Unrecognised varying parameter. Setting to structural discrepancy (disc).")
    modified = 'disc'
  }
  if (modified == 'disc') {
    em_exps <- do.call('cbind', purrr::map(emulators, ~.$get_exp(on_grid)))
    em_vars <- do.call('cbind', purrr::map(emulators, ~.$get_cov(on_grid)))
    #valid_exps <- data.frame(purrr::map(emulators, ~.$get_exp(validation_points[,in_names])))
    #valid_vars <- data.frame(purrr::map(emulators, ~.$get_cov(validation_points[,in_names])))
    #misclass_arr <- array(0, dim = c(length(intervals), length(u_mod)))
    for (i in u_mod) {
      imps <- abs(sweep(em_exps, 2, z_vals, "-"))/sqrt(sweep(em_vars, 2, (i*z_sigs)^2, "+"))
      m_imps <- apply(imps, 1, max)
      cutoff <- purrr::map_dbl(intervals, ~1-length(m_imps[m_imps <= .])/length(m_imps))
      #em_imps <- abs(valid_exps - z_vals)/sqrt(valid_vars + (z_sigs * i)^2)
      #sim_imps <- abs(validation_points[,!names(validation_points) %in% in_names] - z_vals)/(z_sigs * i)
      #misc <- purrr::map_dbl(intervals, ~sum(apply(em_imps > . & sim_imps <= ., 1, purrr::some, isTRUE), na.rm = TRUE))/length(validation_points[,1])
      imp_array[, match(i, u_mod)] <- cutoff
      #misclass_arr[, match(i, u_mod)] <- misc
    }
  }
  else {
    for (i in u_mod) {
      if (modified == 'var')
        ems <- purrr::map(emulators, ~.$set_sigma(i*.$u_sigma))
      else {
        ems <- purrr::map(emulators, ~.$set_theta(i*.$theta))
      }
      imps <- nth_implausible(ems, on_grid, z)
      cutoff <- purrr::map_dbl(intervals, ~1-length(imps[imps <= .])/length(imps))
      imp_array[, match(i, u_mod)] <- cutoff
    }
  }
  df1 <- setNames(data.frame(imp_array), u_mod)
  df1$cutoff <- intervals
  #df2 <- setNames(data.frame(misclass_arr), u_mod)
  #df2$cutoff <- intervals
  melted_df1 <- reshape2::melt(df1, id.vars = 'cutoff')
  #melted_df2 <- reshape2::melt(df2, id.vars = 'cutoff')
  tit <- switch(modified, 'disc' = 'structural discrepancy', 'var' = 'variance inflation', 'corr' = 'correlation length inflation')
  subtit <- switch(modified, 'disc' = '% Structural\nDiscrepancy', 'var' = '% Variance\nInflation', 'corr' = '% Theta\nInflation')
  g <- ggplot(data = melted_df1, aes(x = cutoff, y = value, group = variable, colour = variable)) +
    geom_line(lwd = 1.5) +
    #geom_line(data = plyr::mutate(melted_df2, value = value/max(value)), lwd = 0.5) +
    scale_colour_viridis(discrete = TRUE, option = 'cividis', labels = function(b) {paste0(round(as.numeric(b)*100, 0), "%")}) +
    scale_x_continuous("Implausibility cut-off", labels = function(b) {round(b, 1)}) +
    scale_y_continuous("Removed", labels = function(b) {
      paste0(round(b*100,0),"%")
    } #,
    #sec.axis = sec_axis(~.*max(melted_df2$value), name = "Misclassified", labels = function(b) {
    #  paste0(round(b*100,0), "%")
    #})
    ) +
    labs(title = paste("Space removed as a function of implausibility cut-off and", tit), colour = subtit, x = "Cut-off", y = "% removed") +
    theme_minimal()
  return(g)
  #return(df1)
  #return(list(reduced = df1, misclassed = df2))
}

#' Validation Set Comparisons and Implausibility
#'
#' Creates pairs plots on the set of validation points.
#'
#' Plots are organised as:
#'
#' a) Emulated vs Simulator Output (lower diagonal). The emulator outputs are compared
#' against the simulator outputs. Points whose emulated output lies outside the 3-sigma
#' region of the simulated output are coloured red; those inside are coloured green; a
#' gradient between the two extremes indicates goodness-of-fit;
#'
#' b) Implausibility (upper diagonal). The implausibility for each point is calculated,
#' using the same colour scaling as the lower diagonal.
#'
#' @importFrom GGally ggpairs wrap
#' @importFrom rlang quo_get_expr
#' @import ggplot2
#'
#' @param ems The list of trained emulators
#' @param validation_points The validation set to be plotted
#' @param z The target values for each emulated output
#' @param orig_ranges The original ranges for the input parameters (if desired)
#' @param cb Should a colourblind-friendly palette be used for plots? Default: FALSE
#' @param ... Any additional parameters to be passed to internal functions.
#'
#' @return A data.frame containing the validation points, with goodness-of-fit and implausibility.
#'
#' @export
#'
#' @examples
#' ems <- emulator_from_data(GillespieSIR, c('nS','nI','nR'),
#'  ranges = list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05)),
#'  quadratic = TRUE)
#' targets <- list(
#'  list(val = 281, sigma = 10.43),
#'  list(val = 30, sigma = 11.16),
#'  list(val = 689, sigma = 14.32)
#' )
#' validation_pairs(ems, GillespieValidation, targets)
validation_pairs <- function(ems, validation_points, z, orig_ranges, cb = FALSE, ...) {
  if(missing(orig_ranges))
    orig_ranges <- ems[[1]]$ranges
  em_vals <- data.frame(purrr::map(ems, ~.$get_exp(validation_points[,names(ems[[1]]$ranges)])))
  em_uncert <- data.frame(purrr::map(ems, ~.$get_cov(validation_points[,names(ems[[1]]$ranges)])))
  sim_vals <- validation_points[,!(names(validation_points)%in%names(ems[[1]]$ranges))]
  results <- setNames(cbind(validation_points[,names(ems[[1]]$ranges)], apply(abs(em_vals - sim_vals)/sqrt(em_uncert), 1, max)), c(names(ems[[1]]$ranges), 'bad'))
  results$imps <- nth_implausible(ems, validation_points[,names(ems[[1]]$ranges)], z, ...)
  colourbrks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15)
  colournames <- c(0, '', '', 1, '', '', 2, '', '', 3, '', '', '', 5, '', '', '', 10, 15)
  limfun <- function(data, mapping) {
    ggplot(data = data, mapping = mapping) +
      geom_point(cex = 2) +
      xlim(orig_ranges[[rlang::quo_get_expr(mapping$x)]]) +
      ylim(orig_ranges[[rlang::quo_get_expr(mapping$y)]])
  }
  ifelse(cb, cols <- colourblindcont, cols <- redgreencont)
  g <- ggpairs(results, columns = 1:length(orig_ranges), aes(colour = results[,'bad']), legend = c(1,2),
               title = "Emulator Diagnostics (lower) and Emulator Implausibility (upper)",
               lower = list(continuous = wrap(limfun), mapping = aes(colour = results[,'bad'])),
               upper = list(continuous = wrap(limfun), mapping = aes(colour = results[,'imps'])),
               diag = 'blank', progress = FALSE) +
    scale_colour_gradient2(low = cols$low, mid = cols$mid, high = cols$high, midpoint = 3, breaks = colourbrks, name = "Scale", labels = colournames) +
    theme(legend.position = 'right') +
    theme_minimal()
  print(g)
  return(results)
}

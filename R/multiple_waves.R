#' Plot simulator outputs for multiple waves
#'
#' Plots the simulator results for points at successive waves.
#'
#' The values plotted are the outputs from the simulator; the points passed to it are the
#' points suggested by that wave of emulators. By default, wave 0 is included. A colour
#' scheme is chosen outright for all invocations of this function: it is a 10-colour
#' palette. If more waves are required, then an alternative palette should be selected.
#'
#' @import ggplot2
#'
#' @param wave_points The set of wave points, as a list of data.frames
#' @param z The set of target values for each output
#' @param zero_in Is wave zero included? Default: TRUE
#' @param palette If a larger palette is required, it should be supplied here.
#' @param wave_numbers Which waves to plot. If not supplied, all waves are plotted.
#'
#' @return A ggplot object.
#'
#' @export
#'
simulator_plot <- function(wave_points, z, zero_in = TRUE, palette = NULL, wave_numbers = seq(ifelse(zero_in, 0, 1), length(wave_points)-ifelse(zero_in, 1, 0))) {
  variable <- value <- run <- wave <- val <- sigma <- NULL
  output_names <- names(z)
  sim_runs <- do.call('rbind', purrr::map(wave_numbers, ~data.frame(wave_points[[.+1]][,output_names], wave = .)))
  sim_runs$run <- 1:length(sim_runs[,1])
  melted <- reshape2::melt(sim_runs, id.vars = c('run', 'wave'))
  melted$wave = as.factor(melted$wave)
  if (is.null(palette)) {
    pal <- viridisLite::viridis(10, option = 'plasma', direction = -1)
    pal <- pal[seq_along(pal) %in% (wave_numbers+ifelse(zero_in, 1, 0))]
  }
  else pal <- palette
  obs <- data.frame(variable = names(z), val = purrr::map_dbl(z, ~.$val), sigma = purrr::map_dbl(z, ~.$sigma))
  g <- ggplot(data = melted, aes(x = variable, y = value)) +
    geom_line(aes(group = run, colour = wave)) +
    scale_colour_manual(values = pal) +
    geom_point(data = obs, aes(x = variable, y = val)) +
    geom_errorbar(data = obs, aes(y = val, ymax = val + 3*sigma, ymin = val - 3*sigma), width = 0.1, size = 1.25) +
    labs(title = "Simulator evaluations at wave points") +
    theme_minimal()
  return(g)
}

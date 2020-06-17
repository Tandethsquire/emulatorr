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
  if (is.null(palette)) pal <- viridisLite::viridis(10, option = 'plasma', direction = -1)
  else pal <- palette
  pal <- pal[seq_along(pal) %in% (wave_numbers+ifelse(zero_in, 1, 0))]
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

#' Emulator Variance across waves
#'
#' Plots the emulator variance for each output across emulator waves.
#'
#' It is instructive to look at the change in emulator variance over successive
#' waves, rather than across successive outputs. This function provides a means
#' of doing so quickly for each emulator output.
#'
#' A 2d slice is taken across the input space, where mid-range values of any
#' non-plotted parameters are fixed. The emulator variance (or standard
#' deviation) is then calculated across the two-dimensional subspace for each
#' wave, and for each output.
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_spacer
#' @importFrom purrr %>%
#' @importFrom viridis scale_fill_viridis
#'
#' @param waves A list of lists of \code{\link{Emulator}} objects, corresponding to the waves
#' @param output_names The list of desired outputs to be plotted
#' @param plot_dirs The (two) input parameters to be plotted.
#' @param n_points The number of grid points per plotting dimension. Default: 20
#' @param sd Should the standard deviation be plotted instead of the variance? Default: FALSE
#'
#' @return A list of data.frames, each corresponding to a given output over waves.
#'
#' @export
#'
wave_variance <- function(waves, output_names, plot_dirs = names(waves[[1]][[1]]$ranges)[1:2], n_points = 20, sd = FALSE) {
  if (length(plot_dirs) != 2) stop("Two input directions must be specified.")
  variable <- value <- NULL
  main_ranges <- waves[[1]][[output_names[1]]]$ranges
  grid_ranges <- data.frame(purrr::map(names(main_ranges), function(x) {
    if (x %in% plot_dirs) main_ranges[[x]]
    else rep((main_ranges[[x]][1]+main_ranges[[x]][2])/2, 2)
  })) %>% setNames(names(main_ranges))
  on_grid <- unique(setNames(expand.grid(purrr::map(grid_ranges, ~seq(.[[1]], .[[2]], length.out = n_points))), names(main_ranges)))
  output <- purrr::map(output_names, ~setNames(cbind(on_grid, data.frame(purrr::map(waves, function(x) {
    if (!sd) x[[.]]$get_cov(on_grid)
    else sqrt(x[[.]]$get_cov(on_grid))
  }))), c(names(main_ranges), 1:length(waves))))
  melted_frames <- purrr::map(output, ~reshape2::melt(., id.vars = names(main_ranges))) %>% setNames(output_names)
  plot_list <- list()
  for (i in output_names) {
    dat <- melted_frames[[i]]
    plot_bins <- round(seq(min(dat$value), max(dat$value), length.out = 20))
    g <- ggplot(data = dat, aes(x = dat[,plot_dirs[[1]]], y = dat[,plot_dirs[[2]]], group = variable)) +
      geom_contour_filled(aes(z = value), colour = 'black', breaks = plot_bins) +
      scale_fill_viridis(discrete = TRUE, option = 'plasma', labels = plot_bins, name = ifelse(sd, 'SD[f(x)]', 'Var[f(x)]')) +
      facet_wrap(. ~ variable, nrow = length(waves), labeller = labeller(variable = function(x) paste("Wave",x))) +
      labs(title = paste("Variance for output:", i), x = plot_dirs[1], y = plot_dirs[2]) +
      theme_minimal()
    plot_list[[i]] <- g
  }
  while(length(plot_list)%%length(waves) != 0)
    plot_list[[length(plot_list)+1]] <- patchwork::plot_spacer()
  for (i in 0:(ceiling(length(output_names)/length(waves))-1)) {
    plots <- plot_list[(length(waves)*i+1):(length(waves)*i+3)]
    print(patchwork::wrap_plots(plots))
  }
  return(melted_frames)
}

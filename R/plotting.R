#' Plot Emulator Outputs
#'
#' A wrapper for plotting emulator outputs, for two dimensions. If the input space is greater
#' than 2-dimensional, the mid-range values are chosen for any input beyond the first two.
#' If a list of k emulators are given (e.g. those derived from \code{\link{emulator_from_data}}
#' or \code{\link{full_wave}}), then the result is a kx2 grid of plots. If a single emulator is
#' given, then a single plot is returned.
#'
#' Options for plotting variables are passed via the \code{var} parameter: current choices are
#' \code{'exp'} (Expectation), \code{'var'} (Variance), \code{'imp'} (Implausibility), and
#' \code{'maximp'} (n-th maximum Implausibility). If either of the implausibilities are desired,
#' the \code{targets} parameter must not be \code{NULL}. Bear in mind that n-th maximum
#' implausibility is only permitted if a list of emulators is provided: the \code{...}
#' parameters that can be passed to this function are for optional parameters that can be
#' passed to the \code{\link{nth_implausible}} function (for example, which level of
#' maximum implausibility is wanted).
#'
#' @import ggplot2
#' @importFrom viridis scale_fill_viridis
#' @importFrom reshape2 melt
#' @importFrom purrr %>%
#'
#' @param em A single \code{\link{Emulator}} object, or a list thereof
#' @param var_name The output to plot: options are described above
#' @param npoints The number of lattice points per input direction
#' @param targets The observations. Required if plotting implausibility
#' @param ... Any optional parameters for \code{\link{nth_implausible}}
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#'  ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#'  targets <- list(
#'   list(val = 281, sigma = 10.43),
#'   list(val = 30, sigma = 11.16),
#'   list(val = 689, sigma = 14.32)
#'  )
#'  outputs <- c('nS','nI','nR')
#'  ems <- emulator_from_data(GillespieSIR, outputs, ranges, deltas=rep(0.1, 3), quadratic = TRUE)
#'  t_ems <- purrr::map(seq_along(ems), ~ems[[.]]$adjust(GillespieSIR, outputs[[.]]))
#'  names(t_ems) <- outputs
#'  emulator_plot(t_ems$nI)
#'  emulator_plot(t_ems, var_name = 'var', npoints = 10)
#'  emulator_plot(t_ems, var_name = 'sd', npoints = 10)
#'  emulator_plot(t_ems, var_name = 'imp', targets = targets, npoints = 10)
emulator_plot <- function(em, var_name = 'exp', npoints = 40, targets = NULL, ...) {
  makeGrid <- function(n_points, ranges) {
    grd <- expand.grid(seq(ranges[[1]][1], ranges[[1]][2], length.out = n_points), seq(ranges[[2]][1], ranges[[2]][2], length.out = n_points))
    if (length(ranges) > 2) {
      for (i in 3:length(ranges)) {
        nm <- names(ranges)[i]
        grd[[nm]] <- (ranges[[i]][1]+ranges[[i]][2])/2
      }
    }
    return(grd)
  }
  # If a single emulator, plot it!
  if (var_name == 'exp') title = "Emulator Expectation"
  if (var_name == 'var') title = "Emulator Variance"
  if (var_name == 'sd') title = "Emulator Standard Deviation"
  if (var_name == 'imp') title = "Emulator Implausibility"
  if (class(em)[1] == "Emulator")
  {
    grid <- makeGrid(npoints, em$ranges)
    if (var_name == "exp") {
      outputs <- em$get_exp(grid)
      cols <- 'magma'
    }
    else if (var_name == 'var' || var_name  == 'sd') {
      outputs <- em$get_cov(grid)
      if (var_name == 'sd') outputs <- sqrt(outputs)
      cols <- 'plasma'
    }
    else if (var_name == 'imp' && !is.null(targets)) {
      outputs <- em$implausibility(grid, targets)
      impbrks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15)
      impnames <- c(0, '', '', 1, '', '', 2, '', '', 3, '', '', '', 5, '', '', '', 10, 15)
      included <- c(purrr::map_lgl(impbrks[2:length(impbrks)], ~any(data_grid[,'value'] < .)), TRUE)
    }
    else stop("Could not plot output.")
    data_grid <- data.frame(cbind(grid, outputs)) %>% setNames(c(names(em$ranges), var_name))
    g <- ggplot(data = data_grid, aes(x = data_grid[,names(em$ranges)[[1]]], y = data_grid[,names(em$ranges)[[2]]]))
    if (var_name == 'exp' || var_name == 'var' || var_name == 'sd')
    {
      if (var_name == 'exp') nm <- 'E[f(x)]'
      else if (var_name == 'var') nm <- 'Var[f(x)]'
      else nm <- 'SD[f(x)]'
      bin_vals <- floor(seq(min(0,min(data_grid[,var_name])), ceiling(max(data_grid[,var_name])), length.out = 15))
      g <- g +
        geom_contour_filled(aes(z = data_grid[,var_name]), bins = 15, colour = 'black') +
        scale_fill_viridis(discrete = TRUE, option = cols, name = nm, labels = bin_vals)
    }
    else if (var_name == "imp") {
      cols <- c('#00FF00', '#18FF00', '#31FF00', '#49FF00', '#62FF00', '#7AFF00', '#93FF00', '#ABFF00', '#C4FF00', '#DDFF00',
                '#E0E200', '#E4C600', '#E8AA00', '#EC8D00', '#EF7100', '#F35500', '#F73800', '#FB1C00', '#FF0000')
      g <- g +
        geom_contour_filled(aes(z = data_grid[,var_name]), breaks = impbrks[included], colour = 'black') +
        scale_fill_manual(values = cols[included], labels = impnames[included], guide = guide_legend(reverse = TRUE)) +
        labs(fill = "I")
    }
    return(g +
             labs(title = title) +
             xlab(names(em$ranges)[[1]]) +
             ylab(names(em$ranges)[[2]]) +
             scale_x_continuous(expand = c(0,0)) +
             scale_y_continuous(expand = c(0,0)) +
             theme_minimal())
  }
  # If multiple emulators, want to generate a melted data.frame.
  if (class(em)[1] == 'list' && class(em[[1]])[1] == 'Emulator') {
    grid <- makeGrid(npoints, em[[1]]$ranges)
    data_grid <- grid
    if (var_name == 'maximp' && !is.null(targets)) {
      data_grid <- cbind(data_grid, nth_implausible(em, grid, targets, ...)) %>% setNames(c(names(em[[1]]$ranges), "I"))
      impbrks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15)
      impnames <- c(0, '', '', 1, '', '', 2, '', '', 3, '', '', '', 5, '', '', '', 10, 15)
      cols <- c('#00FF00', '#18FF00', '#31FF00', '#49FF00', '#62FF00', '#7AFF00', '#93FF00', '#ABFF00', '#C4FF00', '#DDFF00',
                '#E0E200', '#E4C600', '#E8AA00', '#EC8D00', '#EF7100', '#F35500', '#F73800', '#FB1C00', '#FF0000')
      included <- c(purrr::map_lgl(impbrks[2:length(impbrks)], ~any(data_grid$I < .)), TRUE)
      g <- ggplot(data = data_grid, aes(x = data_grid[,names(em[[1]]$ranges)[[1]]], y = data_grid[,names(em[[1]]$ranges)[[2]]])) +
        geom_contour_filled(aes(z = I), breaks = impbrks[included], colour = 'black') +
        scale_fill_manual(values = cols[included], labels = impnames[included], guide = guide_legend(reverse = TRUE)) +
        labs(title = "Maximum Implausibility", fill = "I") +
        xlab(names(em[[1]]$ranges)[[1]]) +
        ylab(names(em[[1]]$ranges)[[2]]) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme_minimal()
      return(g)
    }
    for (i in 1:length(em))
    {
      if (var_name == 'exp')
        data_grid <- cbind(data_grid, em[[i]]$get_exp(grid))
      else if (var_name == 'var')
        data_grid <- cbind(data_grid, em[[i]]$get_cov(grid))
      else if (var_name == 'sd')
        data_grid <- cbind(data_grid, sqrt(em[[i]]$get_cov(grid)))
      else if (var_name == 'imp' && !is.null(targets))
        data_grid <- cbind(data_grid, em[[i]]$implausibility(grid, targets[[i]]))
      else stop("Could not plot output.")
    }
    data_grid <- data_grid %>% setNames(c(names(em[[1]]$ranges), names(em))) %>% melt(id.vars = names(em[[1]]$ranges))
    g <- ggplot(data = data_grid, aes(x = data_grid[,names(em[[1]]$ranges)[[1]]], y = data_grid[,names(em[[1]]$ranges)[[2]]]))
    if (var_name == 'exp' || var_name == 'var' || var_name == 'sd') {
      max_val <- max(data_grid$value)
      min_val <- min(data_grid$value)
      bin_vals <- floor(seq(min(0,min(data_grid$value)), ceiling(max(data_grid$value)), length.out = 15))
      g <- g + geom_contour_filled(aes(z = data_grid[,'value']), breaks = bin_vals, colour = 'black')
      ifelse(var_name == 'exp', col <- 'magma', col <- 'plasma')
      if (var_name == 'exp') nm <- 'E[f(x)]' else if (var_name == 'var') nm <- 'Var[f(x)]' else nm <- 'SD[f(x)]'
      g <- g + scale_fill_viridis(discrete = TRUE, option = col, name = nm, labels = bin_vals)
    }
    else {
      impbrks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15)
      impnames <- c(0, '', '', 1, '', '', 2, '', '', 3, '', '', '', 5, '', '', '', 10, 15)
      cols <- c('#00FF00', '#18FF00', '#31FF00', '#49FF00', '#62FF00', '#7AFF00', '#93FF00', '#ABFF00', '#C4FF00', '#DDFF00',
                '#E0E200', '#E4C600', '#E8AA00', '#EC8D00', '#EF7100', '#F35500', '#F73800', '#FB1C00', '#FF0000')
      included <- c(purrr::map_lgl(impbrks[2:length(impbrks)], ~any(data_grid[,'value'] < .)), TRUE)
      g <- g +
        geom_contour_filled(aes(z = data_grid[,'value']), breaks = impbrks[included], colour = 'black') +
        scale_fill_manual(values = cols[included], labels = impnames[included], guide = guide_legend(reverse = TRUE)) +
        labs(fill = "I")
    }
    return(g +
             facet_wrap(~variable, ncol = max(3, floor(sqrt(length(em)))+1)) +
             labs(title = title) +
             xlab(names(em[[1]]$ranges)[[1]]) +
             ylab(names(em[[1]]$ranges)[[2]]) +
             scale_x_continuous(expand = c(0,0)) +
             scale_y_continuous(expand = c(0,0)) +
             theme_minimal())
  }
}

#' Emulator Plots with Outputs
#'
#' Plots emulator outputs across a set of points, with the corresponding observations
#' overlaid (with the appropriate uncertainty).
#' If no points are provided, then \code{npoints} points are uniformly sampled from the
#' input region. Else the provided points are used - for example, if a non-implausible
#' space is known.
#'
#' @importFrom ggplot2 ggplot aes labs geom_line geom_point geom_errorbar
#' @importFrom reshape2 melt
#' @importFrom purrr %>%
#'
#' @param emulators A list of \code{\link{Emulator}} objects
#' @param targets A named list of observations, given in the usual form
#' @param points A list of points at which the emulators should be evaluated. Default: \code{NULL}
#' @param npoints If no points provided, how many input points to evaluate? Default: 1000
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#'  ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#'  outputs <- c('nS','nI','nR')
#'  targets <- list(
#'   nS = list(val = 281, sigma = 10.43),
#'   nI = list(val = 30, sigma = 11.16),
#'   nR = list(val = 689, sigma = 14.32)
#'  )
#'  ems <- emulator_from_data(GillespieSIR, outputs, ranges,
#'   deltas = rep(0.1, 3), quadratic = TRUE)
#'  t_ems <- purrr::map(seq_along(ems), ~ems[[.]]$adjust(GillespieSIR, outputs[.]))
#'  output_plot(t_ems, targets)
output_plot <- function(emulators, targets, points = NULL, npoints = 1000) {
  if (is.null(names(targets))) names(targets) <- names(emulators)
  if (is.null(points)) {
    ranges <- emulators[[1]]$ranges
    points <- data.frame(purrr::map(ranges, ~runif(npoints, .x[1], .x[2])))
  }
  em_exp <- data.frame(purrr::map(emulators, ~.x$get_exp(points))) %>% setNames(names(targets))
  em_exp$run <- 1:length(points[,1])
  em_exp <- reshape2::melt(em_exp, id.vars = 'run')
  target_df <- data.frame(label = names(targets), mu = purrr::map_dbl(targets, ~.x$val), sigma = purrr::map_dbl(targets, ~.x$sigma))
  variable <- value <- run <- mu <- sigma <- label <- NULL
  g <- ggplot(data = em_exp, aes(x = variable, y = value)) +
    geom_line(colour = "purple", aes(group = run), lwd = 1) +
    geom_point(data = target_df, aes(x = label, y = mu), size = 2) +
    geom_errorbar(data = target_df, aes(x = label, y = mu,
                                        ymin = mu-3*sigma,
                                        ymax = mu+3*sigma), width = .1, size = 1.25) +
    labs(title = "Emulator runs vs. Observations")
  return(g)
}

#' Implausibility Plots
#'
#' Generates implausibility plots: either optical depth, or minimum implausibility.
#' Simply a plot wrapper to make later things easier for lattice plots, etc.
#' However, it does have value in and of itself.
#'
#' @import ggplot2
#' @importFrom stats loess
#'
#' @param df A data.frame containing the inputs and implausibilities
#' @param names A list of names for the inputs.
#' @param nvars The number of variables in the plotting region: either 1 or 2.
#' @param min Should minimum implausibility be plotted? Default: FALSE
#' @param ticks Are axis ticks required? Default: FALSE
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#'  ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#'  implaus <- GillespieImplausibility[,c(names(ranges), "I")]
#'  targets <- list(
#'   list(val = 281, sigma = 10.43),
#'   list(val = 30, sigma = 11.16),
#'   list(val = 689, sigma = 14.32)
#'  )
#'  op_depth <- optical_depth(targets, ranges, 20, plot_vars = c('aSI','aIR'), imps = implaus)
#'  plot_implausible(op_depth, c('aSI','aIR'), 2, ticks = TRUE)
#'  min_imp <- min_implausibility(targets, ranges, points_per_dim = 20, imps = implaus)
#'  plot_implausible(min_imp, c('aSI','aIR'), 2, min = TRUE, ticks = TRUE)
plot_implausible <- function(df, names, nvars, min = FALSE, ticks = FALSE) {
  Proportion <- Minimum <- NULL
  if (min)
  {
    impbrks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15)
    impnames <- c(0, '', '', 1, '', '', 2, '', '', 3, '', '', '', 5, '', '', '', 10, 15)
    cols <- c('#00FF00', '#18FF00', '#31FF00', '#49FF00', '#62FF00', '#7AFF00', '#93FF00', '#ABFF00', '#C4FF00', '#DDFF00',
              '#E0E200', '#E4C600', '#E8AA00', '#EC8D00', '#EF7100', '#F35500', '#F73800', '#FB1C00', '#FF0000')
    included <- c(purrr::map_lgl(impbrks[2:length(impbrks)], ~any(df$Minimum < .)), TRUE)
    g <- ggplot(data = df, aes(x = df[,names[1]], y = df[,names[2]])) +
      geom_raster(aes(fill = Minimum), interpolate = TRUE) +
      scale_fill_gradient2(low = '#00FF00', mid = '#AAFF00', high = '#FF0000', midpoint = 2, breaks = impbrks, labels = NULL) +
      labs(x = names[1], y = names[2])
  }
  else if (nvars == 1)
  {
    g <- ggplot(data = df, aes(x = df[,names], y = Proportion, group = 1)) +
      ggplot2::geom_smooth(formula = y~x, method = stats::loess, se = FALSE, color = 'black') +
      ggplot2::coord_cartesian(ylim=c(0,1)) +
      labs(x = names, y = "Optical Depth")
  }
  else if (nvars == 2)
  {
    g <- ggplot(data = df, aes(x = df[,names[1]], y = df[,names[2]])) +
      geom_raster(aes(fill = Proportion), interpolate = TRUE) +
      labs(x = names[1], y = names[2])
  }
  g <- g + theme_minimal()
  if (ticks)
    return(g + theme(axis.text.x = element_text(angle = 90)))
  else
    return(g + theme(legend.position = "none", axis.ticks = element_blank(), axis.text = element_blank()))
}

#' Implausibility Plot Lattice
#'
#' Creates a lattice of plots, consisting of optical depths (both univariate and
#' two-variable projections), and minimum implausibility.
#' It is advisable to pass no more than 6 input parameter names to this function,
#' as the graphs will be too small for meaningful analysis if more are provided.
#'
#' @importFrom purrr map
#' @importFrom utils combn
#' @importFrom cowplot ggdraw draw_plot
#'
#' @param implausibilities A list of implausibilities for input points
#' @param targets A list of observations, given in the usual form
#' @param ranges The ranges of each of the input parameters.
#' @param np The number of bins in each part of the plot grids
#' @param ... Any additional inputs (to be passed to \code{\link{optical_depth}})
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#'  imp_data <- GillespieImplausibility
#'  ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#'  targets <- list(
#'   list(val = 281, sigma = 10.43),
#'   list(val = 30, sigma = 11.16),
#'   list(val = 689, sigma = 14.32)
#'  )
#'  plot_lattice(imp_data[,c(names(ranges), "I")], targets, ranges, np = 20)
plot_lattice <- function(implausibilities, targets, ranges, np = 50, ...) {
  nvars <- length(names(ranges))
  name_set <- t(combn(names(ranges), 2))
  coord_set <- t(combn(seq_along(names(ranges)), 2))
  d_od <- apply(name_set, 1, function(x) optical_depth(targets, ranges, np, plot_vars = x, imps = implausibilities, ...))
  s_od <- map(names(ranges), ~optical_depth(targets, ranges, np, plot_vars = .x, imps = implausibilities, ...))
  d_min <- apply(name_set, 1, function(x) min_implausibility(targets, ranges, np, plot_vars = x, imps = implausibilities, ...))
  plots_s <- map(seq_along(names(ranges)), ~list(plt = plot_implausible(s_od[[.x]], names(ranges)[.x], 1), x = .x-1, y = .x-1))
  plots_d <- map(seq_along(name_set[,1]), ~list(plt = plot_implausible(d_od[[.x]], name_set[.x,], 2), x = coord_set[.x,1]-1, y = coord_set[.x,2]-1))
  plots_m <- map(seq_along(name_set[,1]), ~list(plt = plot_implausible(d_min[[.x]], name_set[.x,], 2, min = TRUE), x = coord_set[.x,2]-1, y = coord_set[.x,1]-1))
  tot_plots <- c(plots_s, plots_d, plots_m)
  plot_out <- ggdraw()
  for (i in 1:length(tot_plots)) {
    to_plot <- tot_plots[[i]]
    plot_out <- plot_out + draw_plot(to_plot$plt, x = to_plot$x/nvars, y = to_plot$y/nvars, width = 1/nvars, height = 1/nvars)
  }
  return(plot_out)
}

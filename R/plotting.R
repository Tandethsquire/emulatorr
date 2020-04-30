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
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_gradientn scale_fill_gradient2 labs xlab ylab facet_wrap
#' @importFrom wesanderson wes_palette
#' @importFrom reshape2 melt
#' @importFrom purrr %>%
#'
#' @param em A single \code{\link{Emulator}} object, or a list thereof
#' @param var The output to plot: options are described above
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
#'  emulator_plot(t_ems, var = 'var', npoints = 20)
#'  emulator_plot(t_ems, var = 'imp', targets = targets, npoints = 20)
emulator_plot <- function(em, var = 'exp', npoints = 40, targets = NULL, ...) {
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
  if (var == 'exp') title = "Emulator Expectation"
  if (var == 'var') title = "Emulator Variance"
  if (var == 'imp') title = "Emulator Implausibility"
  if (class(em)[1] == "Emulator")
  {
    grid <- makeGrid(npoints, em$ranges)
    if (var == "exp") {
      outputs <- em$get_exp(grid)
      cols <- 'Royal1'
    }
    else if (var == 'var') {
      outputs <- em$get_cov(grid)
      cols <- 'Royal2'
    }
    else if (var == 'imp' && !is.null(targets)) {
      outputs <- em$implausibility(grid, targets)
      impbrks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15)
      impnames <- c(0, '', '', 1, '', '', '', '', '', 3, '', '', '', 5, '', '', '', 10, 15)
    }
    else stop("Could not plot output.")
    data_grid <- data.frame(cbind(grid, outputs)) %>% setNames(c(names(em$ranges), var))
    g <- ggplot(data = data_grid, aes(x = data_grid[,names(em$ranges)[[1]]], y = data_grid[,names(em$ranges)[[2]]])) +
      geom_raster(aes(fill = data_grid[,var]), interpolate = TRUE)
    if (var == 'exp' || var == 'var')
    {
      g <- g + scale_fill_gradientn(colours = wes_palette(cols, 20, type = 'continuous'), name = ifelse(var == 'exp', 'E[f(x)]', 'Var[f(x)]'))
    }
    else if (var == "imp") {
      g <- g + scale_fill_gradient2(low = '#00FF00', mid = '#DDFF00', high = '#FF0000', midpoint = 3, breaks = impbrks, labels = impnames, name = "I")
    }
    return(g + labs(title = title) + xlab(names(em$ranges)[[1]]) + ylab(names(em$ranges)[[2]]))
  }
  # If multiple emulators, want to generate a melted data.frame.
  if (class(em)[1] == 'list' && class(em[[1]])[1] == 'Emulator') {
    grid <- makeGrid(npoints, em[[1]]$ranges)
    data_grid <- grid
    if (var == 'maximp' && !is.null(targets)) {
      data_grid <- cbind(data_grid, nth_implausible(em, grid, targets, ...)) %>% setNames(c(names(em[[1]]$ranges), "I"))
      impbrks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15)
      impnames <- c(0, '', '', 1, '', '', '', '', '', 3, '', '', '', 5, '', '', '', 10, 15)
      g <- ggplot(data = data_grid, aes(x = data_grid[,names(em[[1]]$ranges)[[1]]], y = data_grid[,names(em[[1]]$ranges)[[2]]])) +
        geom_raster(aes(fill = I), interpolate = TRUE) +
        scale_fill_gradient2(low = '#00FF00', mid = '#DDFF00', high = '#FF0000', midpoint = 3, breaks = impbrks, labels = impnames) +
        labs(title = "Combined Implausibility") +
        xlab(names(em[[1]]$ranges)[[1]]) +
        ylab(names(em[[1]]$ranges)[[2]])
      return(g)
    }
    for (i in 1:length(em))
    {
      if (var == 'exp')
        data_grid <- cbind(data_grid, em[[i]]$get_exp(grid))
      else if (var == 'var')
        data_grid <- cbind(data_grid, em[[i]]$get_cov(grid))
      else if (var == 'imp' && !is.null(targets))
        data_grid <- cbind(data_grid, em[[i]]$implausibility(grid, targets[[i]]))
      else stop("Could not plot output.")
    }
    data_grid <- data_grid %>% setNames(c(names(em[[1]]$ranges), names(em))) %>% melt(id.vars = names(em[[1]]$ranges))
    g <- ggplot(data = data_grid, aes(x = data_grid[,names(em[[1]]$ranges)[[1]]], y = data_grid[,names(em[[1]]$ranges)[[2]]])) +
      geom_raster(aes(fill = data_grid[,'value']), interpolate = TRUE)
    if (var == "exp") g <- g + scale_fill_gradientn(colours = wes_palette('Royal1',20,'continuous'), name = "E[f(x)]")
    else if (var == 'var') g <- g + scale_fill_gradientn(colours = wes_palette('Royal2',20,'continuous'), name = "Var[f(x)]")
    else {
      impbrks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15)
      impnames <- c(0, '', '', 1, '', '', '', '', '', 3, '', '', '', 5, '', '', '', 10, 15)
      g <- g + scale_fill_gradient2(low = '#00FF00', mid = '#DDFF00', high = '#FF0000', midpoint = 3, breaks = impbrks, labels = impnames, name = 'I')
    }
    return(g + facet_wrap(~variable, ncol = 2) + labs(title = title) + xlab(names(em[[1]]$ranges)[[1]]) + ylab(names(em[[1]]$ranges)[[2]]))
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
    geom_line(colour = "purple", aes(group = run)) +
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
#' @importFrom ggplot2 geom_raster geom_smooth ggplot theme_dark theme element_text element_blank coord_cartesian
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
  Minimum <- Proportion <- NULL
  impbrks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3)
  lbls <- c('', '', '', 1, '', '', '2', '', '', 3)
  if (min)
  {
    g <- ggplot(data = df, aes(x = df[,names[1]], y = df[,names[2]])) +
      geom_raster(aes(fill = Minimum), interpolate = TRUE) +
      scale_fill_gradient2(low = '#00FF00', mid = '#AAFF00', high = '#FF0000', midpoint = 2, breaks = impbrks, labels = NULL) +
      labs(x = names[1], y = names[2])
  }
  else if (nvars == 1)
  {
    g <- ggplot(data = df, aes(x = df[,names], y = Proportion, group = 1)) +
      ggplot2::geom_smooth(formula = y~x, method = stats::loess, se = FALSE, color = 'white') +
      ggplot2::coord_cartesian(ylim=c(0,1)) +
      theme_dark() +
      labs(x = names, y = "Optical Depth")
  }
  else if (nvars == 2)
  {
    g <- ggplot(data = df, aes(x = df[,names[1]], y = df[,names[2]])) +
      geom_raster(aes(fill = Proportion), interpolate = TRUE) +
      labs(x = names[1], y = names[2])
  }
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

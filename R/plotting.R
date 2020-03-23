#' Plot Emulator Outputs
#'
#' A wrapper for plotting emulator outputs, for two dimensions.
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
emulator_plot <- function(em, var = 'exp', npoints = 40, targets = NULL, ...) {
  makeGrid <- function(n_points, ranges) {
    return(expand.grid(seq(ranges[[1]][1], ranges[[1]][2], length.out = n_points), seq(ranges[[2]][1], ranges[[2]][2], length.out = n_points)))
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

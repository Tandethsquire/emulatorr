#' Plot Emulator Outputs
#'
#' A wrapper for plotting emulator outputs, for two dimensions.
#'
#' If the input space is greater
#' than 2-dimensional, the mid-range values are chosen for any input beyond the plotted two,
#' unless specific slice values are chosen (see below).
#'
#' If a list of k emulators are given (e.g. those derived from \code{\link{emulator_from_data}}
#' or \code{\link{full_wave}}), then the result is a kx3 grid of plots. If a single emulator is
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
#' The two parameters \code{params} and \code{fixed_vals} determine the 2d slice that is plotted.
#' The argument \code{params} should be a list containing exactly two elements: the names of the
#' two parameters to plot. The argument \code{fixed_vals} can contain any number of remaining
#' parameters and the values to fix them at: this should be a named list of values. If any
#' parameters do not have specified values, their mid-range values are chosen.
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
#' @param cb Should implausibility plots be coloured as colourblind friendly?
#' @param params Which input parameters should be plotted?
#' @param fixed_vals What should the fixed values of unplotted parameters be?
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
#'  emulator_plot(t_ems, npoints = 10, fixed_vals = list(aSR = 0.01))
#'  emulator_plot(t_ems, npoints = 10, params = c('aSI', 'aSR'), fixed_vals = list(aIR = 0.4))
emulator_plot <- function(em, var_name = 'exp', npoints = 40, targets = NULL, cb = FALSE, params = NULL, fixed_vals = NULL, ...) {
  if (class(em)[1] == "Emulator") rs <- em$ranges
  else if (class(em)[1] == "list" && class(em[[1]])[1] == "Emulator") rs <- em[[1]]$ranges
  if (is.null(params) || length(params) != 2 || any(!params %in% names(rs))) p_vals <- c(1,2)
  else p_vals <- which(names(rs) %in% params)
  makeGrid <- function(n_points, ranges) {
    grd <- expand.grid(seq(ranges[[p_vals[1]]][1], ranges[[p_vals[1]]][2], length.out = n_points), seq(ranges[[p_vals[2]]][1], ranges[[p_vals[2]]][2], length.out = n_points))
    grd <- setNames(grd, names(ranges)[p_vals])
    if (!is.null(fixed_vals) && all(names(fixed_vals) %in% names(ranges))) {
      for (i in 1:length(fixed_vals)) grd[[names(fixed_vals)[i]]] <- fixed_vals[[names(fixed_vals)[i]]]
      used_names <- c(names(ranges)[p_vals], names(fixed_vals))
    }
    else used_names <- names(ranges)[p_vals]
    if (length(used_names) < length(ranges))
    {
      unused_nm_ind <- which(!names(ranges) %in% used_names)
      unused_nms <- names(ranges)[unused_nm_ind]
      for (i in 1:length(unused_nms)) {
        grd[[unused_nms[i]]] <- (ranges[[unused_nm_ind[i]]][1] + ranges[[unused_nm_ind[i]]][2])/2
      }
    }
    grd <- grd[,purrr::map_dbl(names(ranges), ~which(names(grd) == .))]
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
      impbrks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15, 1000)
      impnames <- c(0, '', '', 1, '', '', 2, '', '', 3, '', '', '', 5, '', '', '', 10, 15, '')
      included <- c(purrr::map_lgl(impbrks[2:length(impbrks)], ~any(outputs < .)), TRUE)
    }
    else stop("Could not plot output.")
    data_grid <- data.frame(cbind(grid, outputs)) %>% setNames(c(names(em$ranges), var_name))
    g <- ggplot(data = data_grid, aes(x = data_grid[,names(em$ranges)[[p_vals[1]]]], y = data_grid[,names(em$ranges)[[p_vals[2]]]]))
    if (var_name == 'exp' || var_name == 'var' || var_name == 'sd')
    {
      if (var_name == 'exp') nm <- 'E[f(x)]'
      else if (var_name == 'var') nm <- 'Var[f(x)]'
      else nm <- 'SD[f(x)]'
      bin_vals <- unique(floor(seq(min(0,min(data_grid[,var_name])), ceiling(max(data_grid[,var_name])), length.out = 15)))
      g <- g +
        geom_contour_filled(aes(z = data_grid[,var_name]), breaks = bin_vals, colour = 'black') +
        scale_fill_viridis(discrete = TRUE, option = cols, name = nm, labels = bin_vals)
    }
    else if (var_name == "imp") {
      ifelse(cb, cols <- colourblind, cols <- redgreen)
      g <- g +
        geom_contour_filled(aes(z = data_grid[,var_name]), breaks = impbrks[included], colour = 'black') +
        scale_fill_manual(values = cols[included], labels = impnames[included], guide = guide_legend(reverse = TRUE)) +
        labs(fill = "I")
    }
    return(g +
             labs(title = title) +
             xlab(names(em$ranges)[[p_vals[1]]]) +
             ylab(names(em$ranges)[[p_vals[2]]]) +
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
      impbrks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15, 1000)
      impnames <- c(0, '', '', 1, '', '', 2, '', '', 3, '', '', '', 5, '', '', '', 10, 15, '')
      ifelse(cb, cols <- colourblind, cols <- redgreen)
      included <- c(purrr::map_lgl(impbrks[2:length(impbrks)], ~any(data_grid$I < .)), TRUE)
      g <- ggplot(data = data_grid, aes(x = data_grid[,names(em[[1]]$ranges)[[p_vals[1]]]], y = data_grid[,names(em[[1]]$ranges)[[p_vals[2]]]])) +
        geom_contour_filled(aes(z = I), breaks = impbrks[included], colour = 'black') +
        scale_fill_manual(values = cols[included], labels = impnames[included], guide = guide_legend(reverse = TRUE)) +
        labs(title = "Maximum Implausibility", fill = "I") +
        xlab(names(em[[1]]$ranges)[[p_vals[1]]]) +
        ylab(names(em[[1]]$ranges)[[p_vals[2]]]) +
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
    em_names <- purrr::map_chr(em, ~.$output_name)
    data_grid <- data_grid %>% setNames(c(names(em[[1]]$ranges), em_names)) %>% melt(id.vars = names(em[[1]]$ranges))
    g <- ggplot(data = data_grid, aes(x = data_grid[,names(em[[1]]$ranges)[[p_vals[1]]]], y = data_grid[,names(em[[1]]$ranges)[[p_vals[2]]]]))
    if (var_name == 'exp' || var_name == 'var' || var_name == 'sd') {
      max_val <- max(data_grid$value)
      min_val <- min(data_grid$value)
      bin_vals <- unique(floor(seq(min(0,min(data_grid$value)), ceiling(max(data_grid$value)), length.out = 15)))
      g <- g + geom_contour_filled(aes(z = data_grid[,'value']), breaks = bin_vals, colour = 'black')
      ifelse(var_name == 'exp', col <- 'magma', col <- 'plasma')
      if (var_name == 'exp') nm <- 'E[f(x)]' else if (var_name == 'var') nm <- 'Var[f(x)]' else nm <- 'SD[f(x)]'
      g <- g + scale_fill_viridis(discrete = TRUE, option = col, name = nm, labels = bin_vals)
    }
    else {
      impbrks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15, 1000)
      impnames <- c(0, '', '', 1, '', '', 2, '', '', 3, '', '', '', 5, '', '', '', 10, 15, '')
      ifelse(cb, cols <- colourblind, cols <- redgreen)
      included <- c(purrr::map_lgl(impbrks[2:length(impbrks)], ~any(data_grid[,'value'] < .)), TRUE)
      g <- g +
        geom_contour_filled(aes(z = data_grid[,'value']), breaks = impbrks[included], colour = 'black') +
        scale_fill_manual(values = cols[included], labels = impnames[included], guide = guide_legend(reverse = TRUE)) +
        labs(fill = "I")
    }
    return(g +
             facet_wrap(~variable, ncol = max(3, floor(sqrt(length(em)))+1)) +
             labs(title = title) +
             xlab(names(em[[1]]$ranges)[[p_vals[1]]]) +
             ylab(names(em[[1]]$ranges)[[p_vals[2]]]) +
             scale_x_continuous(expand = c(0,0)) +
             scale_y_continuous(expand = c(0,0)) +
             theme_minimal())
  }
}


#' Emulator Plots with Outputs
#'
#' Plots emulator outputs across a set of points, with the corresponding observations
#' overlaid (with the appropriate uncertainty).
#'
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

#' Create a Plot Lattice
#'
#' Plots a number of emulator plots, all projections of the full-dimensional space.
#'
#' The plots to be included are:
#'
#' One dimensional optical depth plots (the proportion of points that are non-implausible)
#'
#' Two dimensional optical depth plots
#'
#' Two dimensional minimum implausibility plots
#'
#' The 1d optical depth plots are situated on the diagonal, the 2d optical depth plots in the
#' upper triangular elements, and the minimum implausibility plots in the lower triangular
#' elements. To evaluate the quantities, a regular grid with \code{ppd} points per dimension
#' is created, and maximum implausibility is calculated over this grid.
#'
#' @importFrom cowplot plot_grid get_legend
#' @importFrom stats xtabs
#' @import ggplot2
#'
#' @param ems The list of emulators
#' @param targets The corresponding list of targets for the emulators
#' @param ppd The number of points to sample per input dimension. The granularity should be
#' carefully considered for large parameter spaces. Default: 20
#' @param cb Should a colourblind-friendly palette be used for implausibility? Default: FALSE
#'
#' @return A ggplot object
#'
#' @export
#'
plot_lattice <- function(ems, targets, ppd = 20, cb = FALSE) {
  get.count <- function(df, nbins, nms) {
    out_df <- setNames(expand.grid(1:nbins, 1:nbins), nms)
    out_df$Freq <- 0
    for (i in 1:nbins) {
      for (j in 1:nbins) {
        out_df[out_df[,nms[1]] == i & out_df[,nms[2]] == j, 'Freq'] <- sum(apply(df, 1, function(x) x[nms[1]] == i && x[nms[2]] == j))
      }
    }
    x_form <- as.formula(paste("Freq ~ ", nms[1], " + ", nms[2], sep = ""))
    return(xtabs(x_form, data = out_df))
  }
  ranges <- ems[[1]]$ranges
  dim_unif <- purrr::map(ranges, ~seq(.[[1]], .[[2]], length.out = ppd))
  grd <- expand.grid(dim_unif)
  imps <- nth_implausible(ems, grd, targets)
  df <- setNames(data.frame(cbind(grd, imps)), c(names(ranges), "I"))
  bins <- ppd-1
  df_filtered <- subset(df, I <= 3)
  onedimops <- data.frame(do.call('cbind', purrr::map(seq_along(ranges), function(x) {
    binpoints <- seq(ranges[[x]][1], ranges[[x]][2], length.out = bins)
    intervals <- findInterval(df[,names(ranges)[[x]]], binpoints)
    tot_in_bin <- c(tabulate(intervals, nbins = bins))
    tot_NROY <- c(tabulate(findInterval(df_filtered[,names(ranges)[[x]]], binpoints), nbins = bins))
    props <- tot_NROY/tot_in_bin
    return(props[intervals])
  }))) %>% setNames(names(ranges))
  variable_combs <- unlist(purrr::map(outer(names(ranges), names(ranges), paste)[upper.tri(outer(names(ranges), names(ranges), paste))],
                                      ~strsplit(., " ")), recursive = FALSE)
  twodimops <- purrr::map(variable_combs, function(x) {
    binpoints <- purrr::map(seq_along(x), ~seq(ranges[[x[.]]][1], ranges[[x[.]]][2], length.out = bins))
    intervals <- data.frame(do.call('cbind', purrr::map(seq_along(x), ~findInterval(df[,x[.]], binpoints[[.]])))) %>% setNames(x)
    tot_in_bin <- table(intervals)
    NROY_class <- data.frame(do.call('cbind', purrr::map(seq_along(x), ~findInterval(df_filtered[,x[.]], binpoints[[.]])))) %>% setNames(x)
    tot_NROY <- get.count(NROY_class, bins, x)
    props <- tot_NROY/tot_in_bin
    return(apply(intervals, 1, function(x) props[x[1], x[2]]))
  })
  twodimmins <- purrr::map(variable_combs, function(x) {
    binpoints <- purrr::map(seq_along(x), ~seq(ranges[[x[.]]][1], ranges[[x[.]]][2], length.out = bins))
    intervals <- data.frame(do.call('cbind', purrr::map(seq_along(x), ~findInterval(df[,x[.]], binpoints[[.]])))) %>% setNames(x)
    minimp <- matrix(rep(0, bins^2), nrow = bins)
    for (i in 1:bins) {
      for (j in 1:bins) {
        minimp[i, j] <- min(df[intervals[,1] == i & intervals[,2] == j, 'I'])
      }
    }
    return(apply(intervals, 1, function(x) minimp[x[1], x[2]]))
  })
  full_df <- data.frame(cbind(onedimops, twodimops, twodimmins)) %>% setNames(c(paste0(names(ranges), "op"),
                                                                                paste0(purrr::map_chr(variable_combs, ~paste(., collapse="")), "op"),
                                                                                paste0(purrr::map_chr(variable_combs, ~paste(., collapse="")), "min")))
  full_df <- data.frame(cbind(grd, full_df))
  u_mat <- l_mat <- matrix(rep("", length(ranges)^2), nrow = length(ranges))
  comb_names <- purrr::map(variable_combs, paste, collapse = "")
  u_mat[upper.tri(u_mat)] <- paste0(comb_names, "min")
  l_mat[upper.tri(l_mat)] <- paste0(comb_names, "op")
  l_mat <- t(l_mat)
  name_mat <- matrix(paste(u_mat, l_mat, sep = ""), nrow = length(ranges))
  diag(name_mat) <- paste0(names(ranges), 'op')
  ifelse(cb, cols <- colourblind, cols <- redgreen)
  impbrks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15)
  impnames <- c(0, '', '', 1, '', '', 2, '', '', 3, '', '', '', 5, '', '', '', 10, 15)
  plot_list <- unlist(purrr::map(seq_along(1:nrow(name_mat)), function(x) {
    purrr::map(seq_along(1:ncol(name_mat)), function(y) {
      if (x == y) {
        g <- ggplot(data = full_df, aes(x = full_df[,names(ranges)[x]])) +
          geom_smooth(aes(y = full_df[,name_mat[x,x]]), colour = 'black')
      }
      else {
        g <- ggplot(data = full_df, aes(x = full_df[,names(ranges)[y]], y = full_df[,names(ranges)[x]]))
        if (x < y) {
          g <- g +
            geom_contour_filled(aes(z = full_df[,name_mat[x,y]]), breaks = impbrks, colour = 'black') +
            scale_fill_manual(values = cols, labels = impnames, name = "Min. I", guide = guide_legend(reverse = TRUE))
        }
        else {
          g <- g +
            geom_raster(aes(fill = full_df[,name_mat[x,y]]), interpolate = TRUE) +
            scale_fill_gradient(low = 'black', high = 'white', breaks = seq(0, 1, by = 0.1), name = "Op. Depth")
        }
      }
      xlab <- ylab <- NULL
      if (x == nrow(name_mat)) xlab <- names(ranges)[y]
      if (y == 1) ylab <- names(ranges)[x]
      return(g + labs(x = xlab, y = ylab) + theme_minimal())
    })
  }), recursive = FALSE)
  legend_list <- list(get_legend(plot_list[[3]]), rep(NULL, length(ranges)-2), get_legend(plot_list[[10]]))
  for (i in 1:length(plot_list)) plot_list[[i]] <- plot_list[[i]] + theme(legend.position = 'none')
  return(suppressMessages(plot_grid(plot_grid(plotlist = plot_list), plot_grid(plotlist = legend_list, ncol = 1), rel_widths = c(1,0.1))))
}

#' Plot Points of Waves
#'
#' Creates a set of pairs plots of the input points for multiple waves.
#'
#' For each pair of inputs, the points from a succession of waves are plotted,
#' coloured according to the wave number. Histograms of the points for each input
#' (broadly interpreted as marginal distributions for the inputs) are provided
#' on the main diagonal.
#'
#' @importFrom purrr %>%
#' @importFrom GGally ggpairs
#'
#' @param pts_list A list object, whose elements are data.frames of points
#' @param in_names The input dimension names.
#'
#' @return A ggplot object.
#' @export
wave_points <- function(pts_list, in_names, surround = FALSE) {
  wave <- NULL
  out_lst <- list()
  for (i in 0:(length(pts_list)-1))
  {
    out_lst[[i+1]] <- cbind(pts_list[[i+1]][,in_names], rep(i, nrow(pts_list[[i+1]]))) %>% setNames(c(in_names, 'wave'))
  }
  tot_dat <- do.call('rbind', out_lst)
  tot_dat$wave <- as.factor(tot_dat$wave)
  wrapfun <- function(data, mapping) {
    g <- ggplot(data = data, mapping = mapping) +
      geom_point(cex = 1.5)
    if (surround)
      g <- geom_point(cex = 1.5, pch = 1, colour = 'black')
  }
  pal <- viridis::viridis(length(pts_list), option = 'D', direction = -1)
  ggpairs(tot_dat, columns = 1:length(in_names), aes(colour = wave),
          lower = list(continuous = wrapfun),
          upper = 'blank', title = 'Wave points pairs plot', progress = FALSE) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = alpha(pal, 0.5)) +
    theme_bw()
}

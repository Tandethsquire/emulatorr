#' Optical Depth Of Emulator Implausibility
#'
#' Calculates optical depth of the implausibility of a set of emulators, given target values.
#' Implausibility is calculated for each (univariate) emulator and the corresponding target
#' value, and the n-th implausibility is computed. This occurs in a lattice hypercube of size
#' points_per_dim^n, and all points are put into bins in the two target directions. The optical
#' depth is the proportion of points in a given bin that are non-implausible.
#'
#' One or two target directions can be specified.
#'
#' @param targets The target values of the outputs, in the normal form
#' @param ranges The ranges of the input variables, as a named list
#' @param points_per_dim The number of lattice points in each input direction
#' @param plot_vars A list of two names, or a single name, corresponding to the plotting inputs
#' @param emulators A set of \code{\link{Emulator}} objects
#' @param imps If implausibility has already been calculated across a grid, the \code{data.frame} goes here
#' @param cutoff The implausibility cutoff. Default: 3
#' @param ... Any options to be passed to the \code{\link{nth_implausible}} function.
#'
#' @return A \code{data.frame} corresponding to the bins and proportions.
#' @export
#'
#' @examples
#'  ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#'  out_vars <- c('nS', 'nI', 'nR')
#'  targets <- list(
#'   list(val = 281, sigma = 37.26),
#'   list(val = 30, sigma = 11.16),
#'   list(val = 689, sigma = 31.72))
#'  first_ems <- emulator_from_data(GillespieSIR, output_names = out_vars,
#'   ranges = ranges, quadratic = TRUE)
#'  trained_ems <- purrr::map(seq_along(first_ems),
#'   ~first_ems[[.x]]$adjust(GillespieSIR, out_vars[[.x]]))
#'  optical_depth(targets, ranges, 5, emulators = trained_ems)
#'  optical_depth(targets, ranges, 5, plot_vars = c('aIR', 'aSR'),
#'  emulators = trained_ems, cutoff = 4, n = 2, max_imp = 30)
#'  optical_depth(targets, ranges, 5, plot_vars = 'aSR',
#'  emulators = trained_ems, cutoff = 2, n = 2)
optical_depth <- function(targets, ranges, points_per_dim, plot_vars = names(ranges)[1:2], emulators = NULL, imps = NULL, cutoff = 3, ...) {
  dim_seqs <- purrr::map(ranges, ~seq(.x[[1]], .x[[2]], length.out = points_per_dim))
  if (is.null(imps)) eval_grid <- expand.grid(dim_seqs)
  else eval_grid <- imps[,names(ranges)]
  if (is.null(imps)) imp_list <- nth_implausible(emulators, eval_grid, targets, ...)
  else imp_list <- imps$I
  if (length(plot_vars) == 2) {
    centers_1 <- round(seq((dim_seqs[[plot_vars[[1]]]][[1]]+dim_seqs[[plot_vars[[1]]]][[2]])/2, (dim_seqs[[plot_vars[[1]]]][[points_per_dim]]+dim_seqs[[plot_vars[[1]]]][[points_per_dim-1]])/2, length.out = points_per_dim-1),4)
    centers_2 <- round(seq((dim_seqs[[plot_vars[[2]]]][[1]]+dim_seqs[[plot_vars[[2]]]][[2]])/2, (dim_seqs[[plot_vars[[2]]]][[points_per_dim]]+dim_seqs[[plot_vars[[2]]]][[points_per_dim-1]])/2, length.out = points_per_dim-1),4)
    lvls_1 <- cut(eval_grid[[plot_vars[[1]]]], points_per_dim-1, centers_1)
    lvls_2 <- cut(eval_grid[[plot_vars[[2]]]], points_per_dim-1, centers_2)
    full_grid <- setNames(cbind(eval_grid, lvls_1, lvls_2, imp_list), c(names(ranges), paste0(plot_vars, "bin"), "I"))
    unfiltered_table <- table(full_grid[,paste0(plot_vars, "bin")])
    filtered_grid <- full_grid[full_grid$I<=cutoff,]
    filtered_table <- table(filtered_grid[,paste0(plot_vars, "bin")])
    depth <- filtered_table/unfiltered_table
  }
  else {
    centers <- round(seq((dim_seqs[[plot_vars]][[1]]+dim_seqs[[plot_vars]][[2]])/2, (dim_seqs[[plot_vars]][[points_per_dim]]+dim_seqs[[plot_vars]][[points_per_dim-1]])/2, length.out = points_per_dim-1), 4)
    lvls <- cut(eval_grid[[plot_vars]], points_per_dim-1, centers)
    full_grid <- setNames(cbind(eval_grid, lvls, imp_list), c(names(ranges), paste(plot_vars, "bin"), "I"))
    unfiltered_tab <- table(full_grid[,paste(plot_vars,"bin")])
    filtered_grid <- full_grid[full_grid$I<=cutoff,]
    filtered_tab <- table(filtered_grid[,paste(plot_vars,"bin")])
    depth <- filtered_tab/unfiltered_tab
  }
  return(setNames(data.frame(depth), c(plot_vars, 'Proportion')))
}

#' Minimum Implausibilities on a Lattice
#'
#' Calculates minimum implausibility of a set of emulators, given target values.
#' Implausibility is calculated for each (univariate) emulator and the corresponding target
#' value, and the n-th implausibility is computed. This occurs in a lattice hypercube of size
#' points_per_dim^n, and all points are put into bins in the two target directions. The
#' minimum is the smallest such value that occurs in each bin.
#'
#' One or two target directions can be specified.
#'
#' @param targets The target values of the outputs, in the normal form
#' @param ranges The ranges of the input variables, as a named list
#' @param points_per_dim The number of lattice points in each input direction
#' @param plot_vars A list of two names, or a single name, corresponding to the plotting inputs
#' @param emulators A set of \code{\link{Emulator}} objects
#' @param imps If implausibility has already been calculated across a grid, the \code{data.frame} goes here
#' @param ... Any options to be passed to the \code{\link{nth_implausible}} function.
#'
#' @importFrom stats aggregate
#'
#' @return A \code{data.frame} corresponding to the bins and proportions.
#' @export
#'
#' @examples
#'  ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#'  out_vars <- c('nS', 'nI', 'nR')
#'  targets <- list(
#'   list(val = 281, sigma = 37.26),
#'   list(val = 30, sigma = 11.16),
#'   list(val = 689, sigma = 31.72))
#'  first_ems <- emulator_from_data(GillespieSIR, output_names = out_vars,
#'   ranges = ranges, quadratic = TRUE)
#'  trained_ems <- purrr::map(seq_along(first_ems),
#'   ~first_ems[[.x]]$adjust(GillespieSIR, out_vars[[.x]]))
#'  min_implausibility(targets, ranges, 5, emulators = trained_ems)
#'  min_implausibility(targets, ranges, 5, plot_vars = c('aIR', 'aSR'),
#'  emulators = trained_ems, n = 2, max_imp = 30)
min_implausibility <- function(targets, ranges, points_per_dim, plot_vars=names(ranges)[1:2], emulators = NULL, imps = NULL, ...) {
  dim_seqs <- purrr::map(ranges, ~seq(.x[[1]], .x[[2]], length.out = points_per_dim))
  if (is.null(imps)) eval_grid <- expand.grid(dim_seqs)
  else eval_grid <- imps[,names(ranges)]
  if (is.null(imps)) imp_list <- nth_implausible(emulators, eval_grid, targets, ...)
  else imp_list <- imps$I
  centers_1 <- round(seq((dim_seqs[[plot_vars[[1]]]][[1]]+dim_seqs[[plot_vars[[1]]]][[2]])/2, (dim_seqs[[plot_vars[[1]]]][[points_per_dim]]+dim_seqs[[plot_vars[[1]]]][[points_per_dim-1]])/2, length.out = points_per_dim-1),4)
  centers_2 <- round(seq((dim_seqs[[plot_vars[[2]]]][[1]]+dim_seqs[[plot_vars[[2]]]][[2]])/2, (dim_seqs[[plot_vars[[2]]]][[points_per_dim]]+dim_seqs[[plot_vars[[2]]]][[points_per_dim-1]])/2, length.out = points_per_dim-1),4)
  lvls_1 <- cut(eval_grid[[plot_vars[[1]]]], points_per_dim-1, centers_1)
  lvls_2 <- cut(eval_grid[[plot_vars[[2]]]], points_per_dim-1, centers_2)
  full_grid <- setNames(cbind(eval_grid, lvls_1, lvls_2, imp_list), c(names(ranges), paste0(plot_vars, 'bin'), 'I'))
  mins <- aggregate(as.formula(paste("I~",paste0(paste(plot_vars,'bin',sep=""),collapse=" + "),sep=" ")), data = full_grid, function(x) min(x))
  return(setNames(mins, c(plot_vars, 'Minimum')))
}

#' Generate Simulator Runs
#'
#' A wrapper for a variety of sampling methods.
#' Given a set of trained emulators, finds the next set of points that will be informative
#' for the next wave of emulators.
#'
#' If the method is 'lhs', this creates a new training set using LHS, and
#' then finds the trace of the variance matrix of the emulators across these points (this is
#' broadly equivalent to using V-optimality). We repeat this for \code{n_runs}, and select
#' the configuration that minimises the mean of the variances across the emulators.
#' If observations are given, then these are used to ensure that the new sample points are
#' non-implausible according to the current emulators.
#'
#' If the method is 'slice', then an initial point \code{x} must be provided. This is ideally
#' an non-implausible point, or at least one close to the suspected non-implausible region.
#' It then applies slice sampling, using implausibility as a measure of success.
#'
#' If the method is 'optical', then the optical depth of the space in each parameter direction
#' is calculated (using a known set of non-implausible points \code{plausible_set}),
#' and used as a distribution for that parameter. Points are sampled from the
#' collection of distributions and non-implausible points generated are filtered out. From the
#' remaining points, a sample of the required size is generated using maximin criterion.
#'
#' For any sampling strategy, the parameters \code{emulators} and \code{ranges} must be
#' specified. If the method is 'slice', then the parameters \code{x}
#' and \code{z} are necessary. All other parameters are optional.
#'
#' @importFrom stats setNames runif dist
#'
#' @param emulators A list of \code{\link{Emulator}} objects, trained on the design points
#' @param ranges The ranges of the input parameters
#' @param n_points Optional. Specifies how many additional points are required. Default: 10*(number of emulators)
#' @param z Checks implausibility of sample points to restrict to only non-implausible points.
#' @param method Any of 'lhs', 'slice', 'optical'.
#' @param cutoff Optional. If z is given, this is the implausibility cutoff for the filtering. Default = 3
#' @param ... Any parameters that need to be passed to a particular method (see above)
#'
#' @return A \code{data.frame} containing the set of new points to simulate at.
#' @export
#'
#' @examples
#' ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#' ems <- emulator_from_data(GillespieSIR, output_names = c('nS', 'nI', 'nR'),
#'  ranges = ranges, quadratic = TRUE)
#' trained_ems <- purrr::map(seq_along(ems),
#'  ~ems[[.x]]$adjust(GillespieSIR, c('nS', 'nI', 'nR')[[.x]]))
#' targets <- list(
#'  list(val = 281, sigma = 10.43),
#'  list(val = 30, sigma = 11.16),
#'  list(val = 689, sigma = 14.32)
#' )
#' non_imp_points <- GillespieImplausibility[GillespieImplausibility$I2 <= 4, names(ranges)]
#' example_point <- unlist(non_imp_points[sample(1:length(non_imp_points[,1]), 1),],
#'  use.names = FALSE)
#' pts_lhs <- generate_new_runs(trained_ems, ranges, 10, targets, cutoff = 3)
#' pts_slice <- generate_new_runs(trained_ems, ranges, 10, targets,
#'  method = 'slice', cutoff = 4, x = example_point)
#' pts_optical <- generate_new_runs(trained_ems, ranges, 10, targets,
#'  method = 'optical', cutoff = 4, plausible_set = non_imp_points)
generate_new_runs <- function(emulators, ranges, n_points = 10*length(ranges), z, method = 'lhs', cutoff = 3, ...) {
  if (method == 'lhs')
    return(lhs_generation(emulators, ranges, n_points, z, cutoff, ...))
  if (method == 'slice')
    return(slice_generation(emulators, ranges, n_points, z, cutoff, ...))
  if (method == 'optical')
    return(optical_depth_generation(emulators, ranges, n_points, z, cutoff, ...))
  stop("Unknown method used. Options are 'lhs', 'slice', 'optical'.")
}

# A function to perform lhs sampling
lhs_generation <- function(emulators, ranges, n_points, z, n_runs = 20, cutoff = 3) {
  current_trace = NULL
  out_points <- NULL
  new_points <- eval_funcs(scale_input, setNames(data.frame(2*(lhs::maximinLHS(n_points*20, length(ranges))-1/2)), names(ranges)), ranges, FALSE)
  if (!missing(z)) new_points <- new_points[nth_implausible(emulators, new_points, z)<=cutoff,]
  if (length(new_points[,1]) < n_points) {
    message(cat("Only", length(new_points[,1]), "points generated."))
    return(setNames(data.frame(new_points), names(ranges)))
  }
  else message(cat(length(new_points[,1]), "non-implausible points generated. Applying V-optimality..."))
  for (i in 1:n_runs)
  {
    new_points <- new_points[sample(seq_along(new_points[,1]), n_points),]
    measure <- mean(purrr::map_dbl(seq_along(emulators), ~sum(emulators[[.x]]$get_cov(new_points))))
    if (is.null(current_trace) || measure < current_trace) {
      out_points <- new_points
      current_trace <- measure
    }
  }
  return(setNames(data.frame(out_points), names(ranges)))
}

# A function to perform slice sampling
slice_generation <- function(emulators, ranges, n_points, z, cutoff = 3, x) {
  out_points <- x
  x0 <- c(out_points, use.names = F)
  x_new = x0
  for (i in 1:n_points) {
    for (j in 1:length(ranges)) {
      xl = ranges[[j]][1]
      xr = ranges[[j]][2]
      repeat {
        x_new[j] <- runif(1, min = xl, max = xr)
        if (nth_implausible(emulators, x_new, z)<=cutoff) break
        else ifelse(x_new[j]<x0[j], xl <- x_new[j], xr <- x_new[j])
      }
    }
    out_points <- cbind(out_points, x_new)
    x0 <- x_new
  }
  return(setNames(data.frame(t(out_points), row.names = NULL), names(ranges)))
}

# A function to perform optical depth sampling
optical_depth_generation <- function(emulators, ranges, n_points, z, n_runs = 100, cutoff = 3, plausible_set, ...) {
  get_depth <- function(plausible_set, ranges, var_name, nbins = 100) {
    output <- c()
    varseq <- seq(ranges[[var_name]][1], ranges[[var_name]][2], length.out = (nbins+1))
    for (i in 1:(length(varseq)-1)) {
      opdepth <- length(plausible_set[plausible_set[[var_name]]>=varseq[i] & plausible_set[[var_name]]<varseq[i+1],1])/length(plausible_set[,1])
      output <- c(output, (varseq[i]+varseq[i+1])/2, opdepth)
    }
    return(setNames(data.frame(t(matrix(output, nrow=2))), c('bin', 'prob')))
  }
  output <- c()
  for (i in 1:length(ranges)) {
    probs <- get_depth(plausible_set, ranges, names(ranges)[i], ...)
    new_pts <- sample(probs$bin, n_points*10, prob = probs$prob, replace = TRUE) + runif(n_points*10, min = 0, max = probs$bin[2]-probs$bin[1])
    output <- c(output, new_pts)
  }
  df <- setNames(data.frame(matrix(output, nrow = n_points*10)), names(ranges))
  df <- df[nth_implausible(emulators, df, z) <= cutoff,]
  if (length(df[,1]) < n_points)
  {
    message(cat("Only", length(df[,1]), "points generated."))
    return(df)
  }
  cdist <- 0
  for (i in 1:n_runs) {
    tempset <- df[sample(seq_along(df[,1]), n_points),]
    tempdist <- min(c(dist(tempset)))
    if (tempdist > cdist) {
      outset <- tempset
      cdist <- tempdist
    }
  }
  message(cat("Minimum distance:", cdist))
  return(outset)
}

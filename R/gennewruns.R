# Helper Functions for uniform sphere sampling
#'
#' @importFrom stats rexp rnorm
runifs <- function(n, d, c = rep(0, d), r = 1) {
  init_points <- matrix(rnorm(n*d, 0, 1), nrow = n, byrow = T)
  exps <- rexp(n, 0.5)
  denoms <- (exps + apply(init_points, 1, function(x) sum(x^2)))^(0.5)
  swept <- sweep(init_points, 1, denoms, "/")
  return(t(apply(swept, 1, function(x) x*r + c)))
}
punifs <- function(x, c = rep(0, length(x)), r = 1) {
  return(ifelse(sum((x-c)^2) <= r^2, 1, 0))
}


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
#' If the method is 'importance', importance sampling is used. Starting from a set of
#' non-implausible (preferably space-filling) points, points are sampled from a distribution
#' around the points, and included in the output based on a weighted measure gained from the
#' mixture distribution of the initial points. The set \code{plausible_set} must be specified.
#' If \code{burn_in} is TRUE, then a burn-in phase is used to determine the optimal parameters
#' for the proposal distribution.
#'
#' Note that the \code{plausible_set} parameter size differs between the two methods that use
#' it (the 'optical' and 'importance') methods. The optical set should be as large as possible
#' in order to accurately represent the optical depth in each parameter direction; the set for
#' importance sampling should be smaller (and probably smaller than the desired number of
#' output points) in order to expedite the initial set-up of the sampling strategy.
#'
#' For any sampling strategy, the parameters \code{emulators}, \code{ranges} and
#' \code{z} must be specified. If the method is 'slice', then the parameter \code{x}
#' is necessary. All other parameters are optional.
#'
#' After the first round of sampling, if \code{line_sample} is enabled, an exploration of the
#' boundary of the non-implausible region is performed as follows. Two points are chosen at
#' random, and a number of points are sampled uniformly along the line connecting these points.
#' The sampled points are tested for implausibility, and (provided more than 50% but not all
#' of the points are non-implausible) the most separated of the points replace the two initial
#' points.
#'
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom stats setNames runif dist
#'
#' @param emulators A list of \code{\link{Emulator}} objects, trained on the design points
#' @param ranges The ranges of the input parameters
#' @param n_points Optional. Specifies how many additional points are required. Default: 10*(number of emulators)
#' @param z Checks implausibility of sample points to restrict to only non-implausible points.
#' @param method Any of 'lhs', 'slice', 'optical'.
#' @param cutoff Optional. If z is given, this is the implausibility cutoff for the filtering. Default = 3
#' @param include_line Should line sampling be applied after point generation? Default: TRUE.
#' @param plausible_set Optional - a set of non-implausible points from which to start.
#' @param burn_in If importance sampling, should a burn-in phase be used? Default: FALSE
#' @param ... Any parameters that need to be passed to a particular method (see below)
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
#' non_imp_points <- GillespieImplausibility[GillespieImplausibility$I <= 4, names(ranges)]
#' \donttest{
#' pts_lhs <- generate_new_runs(trained_ems, ranges, 10, targets, cutoff = 3)
#' pts_slice <- generate_new_runs(trained_ems, ranges, 10, targets,
#'  method = 'slice', cutoff = 4, plausible_set = non_imp_points, include_line = FALSE)
#' pts_optical <- generate_new_runs(trained_ems, ranges, 10, targets,
#'  method = 'optical', cutoff = 4, plausible_set = non_imp_points, include_line = FALSE)
#' non_imp_sample <- non_imp_points[sample(seq_along(non_imp_points[,1]), 20),]
#' pts_importance <- generate_new_runs(trained_ems, ranges, 10, targets,
#'  method = 'importance', cutoff = 4, plausible_set = non_imp_sample, include_line = FALSE)}
generate_new_runs <- function(emulators, ranges, n_points = 10*length(ranges), z, method = 'importance', include_line = TRUE, cutoff = 3, plausible_set, burn_in = FALSE, ...) {
  if (missing(plausible_set) || method == 'lhs')
    points <- lhs_generation(emulators, ranges, n_points, z, cutoff)
  else points <- plausible_set
  if (method == 'lhs') return(points)
  if (include_line) points <- line_sample(emulators, points, z, ranges)
  if (method == 'importance') {
    if (!burn_in)
      sd <- min(purrr::map_dbl(ranges, ~.[[2]]-.[[1]]))/8
    else
      sd <- NULL
    points <- importance_sample(emulators, ranges, n_points, z, cutoff, sd, plausible_set = points, ...)
  }
  else if (method == 'slice')
    points <- slice_generation(emulators, ranges, n_points, z, cutoff, x = unlist(points[sample(nrow(points), 1),]), ...)
  else if (method == 'optical')
    points <- optical_depth_generation(emulators, ranges, n_points, z, cutoff, plausible_set = points)
  else stop("Unknown method used. Options are 'lhs', 'slice', 'optical', 'importance'.")
  return(points)
}

# A function to perform lhs sampling
lhs_generation <- function(emulators, ranges, n_points, z, n_runs = 20, cutoff = 3) {
  print("Performing Latin hypercube sampling with rejection...")
  current_trace = NULL
  out_points <- NULL
  new_points <- eval_funcs(scale_input, setNames(data.frame(2*(lhs::maximinLHS(n_points*20, length(ranges))-1/2)), names(ranges)), ranges, FALSE)
  if (!missing(z)) new_points <- new_points[nth_implausible(emulators, new_points, z)<=cutoff,]
  if (length(new_points[,1]) < n_points) {
    message(cat("Only", length(new_points[,1]), "points generated."))
    return(setNames(data.frame(new_points), names(ranges)))
  }
  #else message(cat(length(new_points[,1]), "non-implausible points generated. Applying V-optimality..."))
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
  print("Performing slice sampling...")
  J <- function(x) {
    for (i in 1:length(emulators)) {
      if (!emulators[[i]]$implausibility(x, z[[i]], cutoff)) return(FALSE)
    }
    return(TRUE)
  }
  out_points <- x
  x0 <- x_new <- c(out_points, use.names = F)
  for (i in 1:n_points) {
    for (j in 1:length(ranges)) {
      xl = ranges[[j]][1]
      xr = ranges[[j]][2]
      repeat {
        x_new[j] <- runif(1, min = xl, max = xr)
        if (J(data.frame(matrix(x_new, nrow = 1)) %>% setNames(names(ranges)))) break
        ifelse(x_new[j]<x0[j], xl <- x_new[j], xr <- x_new[j])
      }
    }
    out_points <- cbind(out_points, x_new)
    x0 <- x_new
  }
  return(setNames(data.frame(t(out_points), row.names = NULL), names(ranges))[2:n_points+1,])
}

# A function to perform optical depth sampling
optical_depth_generation <- function(emulators, ranges, n_points, z, n_runs = 100, cutoff = 3, plausible_set, ...) {
  print("Performing optical depth sampling...")
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
    message(cat("Only", length(df[,1]), "points generated from LHS with rejection."))
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

# Importance Sampling from initial points
importance_sample <- function(ems, ranges, n_points, z, cutoff = 3, sd = NULL, dist = 'sphere', plausible_set, burn_in = FALSE) {
  print("Performing importance sampling...")
  J <- function(x) {
    for (i in 1:length(ems)) {
      if (!ems[[i]]$implausibility(x, z[[i]], cutoff)) return(FALSE)
    }
    return(TRUE)
  }
  range_func <- function(x, ranges) {
    all(purrr::map_lgl(seq_along(ranges), ~x[.]>=ranges[[.]][1] && x[.]<=ranges[[.]][2]))
  }
  if (!burn_in) sd = min(purrr::map_dbl(ranges, ~.[[2]]-.[[1]]))/8
  if (is.null(sd)) {
    if (dist == 'sphere') {
      cat("Performing burn-in...\n")
      bi <- min(purrr::map_dbl(ranges, ~.[[2]]-.[[1]]))/6
      b_int <- bi/10
      bv1 <- importance_sample(ems, ranges, max(500, floor(n_points/5)), z, sd = bi, plausible_set = plausible_set, burn_in = TRUE)
      bv2 <- bv1
      while (bv1 <= bv2) {
        bi <- bi + b_int
        bv1 <- bv2
        bv2 <- importance_sample(ems, ranges, max(500, floor(n_points/5)), z, sd = bi, plausible_set = plausible_set, burn_in = TRUE)
      }
      sd <- bi
      cat("Optimal radius: ", sd, "; approximate volume: ", bv2, "\n", sep = "")
    }
    else {
      lengths <- purrr::map_dbl(ranges, ~.[2]-.[1])
      sd <- diag((lengths/6)^2, length(lengths))
    }
  }
  pts_added <- 0
  pts_tested <- 0
  d <- length(plausible_set)
  out_arr <- array(0, dim = c(n_points, length(plausible_set)))
  checking_points <- plausible_set
  n_initial <- length(checking_points[,1])
  if (dist == 'normal') {
    weights <- apply(checking_points, 1, function(x) 1/n_initial * sum(apply(checking_points, 1, function(y) dmvnorm(x, mean = y, sigma = sd))))
    min_w <- min(weights)
  }
  while (pts_added < n_points) {
    pts_tested <- pts_tested + 1
    which_point <- sample(1:n_initial, 1)
    if (dist == 'normal')
      new_point <- rmvnorm(1, mean = unlist(checking_points[which_point,], use.names = F), sigma = sd)
    else {
      new_point <- runifs(1, d, unlist(checking_points[which_point,], use.names = F), r = sd)
    }
    new_point <- data.frame(matrix(new_point, nrow = 1)) %>% setNames(names(ranges))
    if (!J(new_point) || !range_func(new_point, ranges)) next
    if (dist == 'normal') {
      point_weight <- 1/n_initial * sum(apply(checking_points, 1, function(x) dmvnorm(new_point, mean = x, sigma = sd)))
      choosing_prob <- min_w/point_weight
    }
    else {
      point_weight <- sum(apply(checking_points, 1, function(x) punifs(new_point, x, r = sd)))
      choosing_prob <- 1/point_weight
    }
    pick <- runif(1)
    if (pick <= choosing_prob) {
      pts_added <- pts_added + 1
      out_arr[pts_added,] <- unlist(new_point)
    }
  }
  if (burn_in) return(nrow(checking_points) * n_points/pts_tested * pi^(d/2) * sd^d / gamma(d/2+1))
  return(setNames(data.frame(out_arr), names(plausible_set)))
}

# Line sampling
line_sample <- function(ems, points, z, ranges, n_lines = 20, ppl = 25, cutoff = 3) {
  print("Performing line sampling...")
  out_pts <- points
  range_func <- function(x, ranges) {
    all(purrr::map_lgl(seq_along(ranges), ~x[.]>=ranges[[.]][1] && x[.]<=ranges[[.]][2]))
  }
  i <- 1
  while (i <= n_lines) {
    s_pts <- points[sample(nrow(points), 2),]
    x1 <- s_pts[1,]
    x2 <- s_pts[2,]
    lambda <- seq(-1,1,length.out = ppl)
    enough_pts <- FALSE
    how_many_attempts <- 0
    while(!enough_pts) {
      how_many_attempts <- how_many_attempts + 1
      line_points <- do.call('rbind', purrr::map(lambda, ~x1+.*(x1-x2)))
      line_imps <- nth_implausible(ems, line_points, z)
      pts_bool <- (line_imps <= cutoff & apply(line_points, 1, range_func, ranges))
      rest_pts <- line_points[pts_bool,]
      if (pts_bool[1] || pts_bool[length(pts_bool)]) lambda <- 1.1*lambda
      else if (nrow(rest_pts) < 0.25*nrow(line_points)) lambda <- 0.9*lambda
      else {
        enough_pts <- TRUE
        add_pts <- c()
        for (j in which(pts_bool)) {
          if (!(pts_bool[j-1] && pts_bool[j+1])) add_pts <- c(add_pts, j)
        }
        out_pts <- rbind(out_pts, line_points[add_pts,])
      }
      if(how_many_attempts > 10) break
    }
    if (how_many_attempts > 10) i <- i - 1
    i <- i + 1
  }
  return(out_pts)
}

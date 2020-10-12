#' Crossover step
#'
#' Takes two points, chooses a pivot along the points, and exchanges
#' the values after the pivot. Concretely, given x = (x1, x2, ..., xd) and
#' y = (y1, y2, ... , yd), a value k in 1:d is chosen, and the resulting points
#' are (x1, x2, ..., xk, y(k+1), ... yd) and (y1, y2, ..., yk, x(k+1), xd).
#' If the resulting points are non-implausible, they are kept.
#'
#' @param x The current set of points
#' @param imp_func An implausibility function (defined in \code{IDEMC})
#' @param imps The list of implausibilities for the ladder.
#'
#' @keywords internal
#' @noRd
#'
#' @return A data.frame containing the resulting points (with the crossover, if it passed)
crossover <- function(x, imp_func, imps) {
  n <- nrow(x)
  d <- ncol(x)
  indices <- sample(n, 2)
  c_point <- sample(d-1, 1)
  yi <- unlist(c(x[indices[1], 1:c_point], x[indices[2], (c_point+1):d], use.names = FALSE))
  yj <- unlist(c(x[indices[2], 1:c_point], x[indices[1], (c_point+1):d], use.names = FALSE))
  if (imp_func(setNames(data.frame(matrix(yi, nrow = 1)), names(x)), imps[indices[1]]) && imp_func(setNames(data.frame(matrix(yj, nrow = 1)), names(x)), imps[indices[2]])) {
    x[indices[1],] <- yi
    x[indices[2],] <- yj
  }
  return(x)
}

#' Exchange step
#'
#' Takes two points from neighbouring rungs of the implausibility
#' ladder and swaps them. The changed ladders are kept if the point now
#' moved into a lower implausibility ladder does indeed satisfy the conditions.
#'
#' @param x The current set of points
#' @param imp_func An implausibility function (defined in \code{IDEMC})
#' @param imps The list of implausibilities for the ladder.
#'
#' @keywords internal
#' @noRd
#'
#' @return A data.frame containing the resulting points (with the exchange, if it passed)
exchange <- function(x, imp_func, imps) {
  n <- nrow(x)
  i <- sample(n, 1)
  if (i !=n && i != 1)
  {
    r <- runif(1)
    if (r < 0.5) {
      j <- i
      i <- j - 1
    }
    else j <- i + 1
  }
  else if (i == n) {
    j = n
    i = n-1
  }
  else j <- 2
  if (imp_func(x[i,], imps[j])) {
    temp <- x[i,]
    x[i,] <- x[j,]
    x[j,] <- temp
  }
  return(x)
}

#' Mutation step
#'
#' Takes a point and mutates it using Metropolis-Hastings. For a given rung of the
#' implausibility ladder, we have a partition of the space using k-means clustering,
#' from which we can find the within-cluster variance structure and the complete
#' variance structure. These are used as part of the proposal distribution for new
#' points, y.
#'
#' @param x The current points
#' @param specs The cluster specifications for the ladder rungs
#' @param imp_func An implausibility function (defined in \code{IDEMC})
#' @param imps The list of implausibilities for the ladder
#' @param w The weighting of in-cluster and whole-space proposal distributions
#' @param index The index of the rung to operate on. Defaults to \code{NULL}
#'
#' @keywords internal
#' @noRd
#'
#' @return A data.frame with the resulting points (including mutation if successful)
mutate <- function(x, specs, imp_func, imps, w = 0.8, index = NULL) {
  if(is.null(index)) index <- sample(length(imps), 1)
  chromosome <- unlist(x[index,], use.names = FALSE)
  spec <- specs[[index]]
  imp <- imps[index]
  if (index == 1) {
    y <- unlist(purrr::map_dbl(specs[[1]]$ranges, ~runif(1, .[1], .[2])))
    x[index, ] <- y
    return(x)
  }
  which_c <- function(x) {
    mahal_dists <- purrr::map_dbl(seq_along(spec$vj), ~(x-spec$mu[.,]) %*% minv(spec$vj[[.]]) %*% (x - spec$mu[.,]))
    return(which.min(mahal_dists))
  }
  y <- c(w*mvtnorm::rmvnorm(1, chromosome, spec$vj[[which_c(chromosome)]]) + (1-w)*mvtnorm::rmvnorm(1, chromosome, spec$whole))
  if (imp_func(setNames(data.frame(matrix(y, nrow = 1)), names(x)), imps[index])) {
    if (which_c(y) == which_c(chromosome)) x[index,] <- y
    else {
      num <- w*mvtnorm::dmvnorm(chromosome, y, specs$vj[[which_c(chromosome)]]) + (1-w)*mvtnorm::dmvnorm(chromosome, y, spec$whole)
      denom <- w*mvtnorm::dmvnorm(y, chromosome, spec$vj[[which_c(chromosome)]]) + (1-w)*mvtnorm::dmvnorm(y, chromosome, spec$whole)
      if(runif(1) < num/denom) x[index, ] <- y
    }
  }
  return(x)
}

#' Clustering function
#'
#' Obtains the clusters for a given set of points. BIC is used to determine the
#' optimal number of clusters; once this has been found, the points are put into
#' their clusters and a within-cluster variance matrix is determined. The full
#' space also has its variance calculated.
#'
#' @importFrom stats cov kmeans
#'
#' @param x The set of points
#' @param max_clusters Determines the largest possible set of clusters
#'
#' @keywords internal
#' @noRd
#'
#' @return A list containing cluster centres, variances, and the whole-space variance.
get_specs <- function(x, max_clusters = 10) {
  # For some reason, this works outside of the package, but not inside it.
  # logLik.kmeans = function(object, ...) structure(
  #   object$tot.withinss,
  #   df = nrow(object$centers)*ncol(object$centers),
  #   nobs = length(object$cluster),
  #   class = 'logLik'
  # )
  kBIC <- function(fit) {
    D <- fit$tot.withinss
    df <- nrow(fit$centers) * ncol(fit$centers)
    nobs <- length(fit$cluster)
    return(D + df * log(nobs))
  }
  # If the commented code works, replace kBIC with BIC.
  bics <- purrr::map_dbl(1:max_clusters, ~kBIC(kmeans(x, ., nstart = 20)))
  cents <- which.min(bics)
  clust <- kmeans(x, centers = cents, nstart = 20)
  vmeans <- clust$centers
  vmats <- purrr::map(1:cents, ~cov(x[clust$cluster == .,]))
  vwhole <- cov(x)
  return(list(mu = vmeans, vj = vmats, whole = vwhole))
}

#' IDEMC step
#'
#' Performs a single IDEMC step: given a set of points in an implauisbility ladder,
#' it performs some number of mutations, crossovers, and exchanges. Only returns if
#' a new point has been found.
#'
#' @param points The initial set of points
#' @param imp_func An implausibility function (defined in \code{IDEMC})
#' @param imp_levels The implausibility ladder levels
#' @param all_specs The cluster specifications for each rung of the ladder
#' @param s The number of points to generate
#' @param pm The probability of doing mutation rather than crossover in the step (default: 0.9)
#' @param M If mutation, how many mutations should be performed (default: 10)
#'
#' @keywords internal
#' @noRd
#'
#' @return A set of laddered points generated by the step.
IDEMC_step <- function(points, imp_func, imp_levels, all_specs, s, pm = 0.9, M = 10) {
  total_set <- purrr::map(seq_along(points[,1]), ~points[.,])
  start_points <- points
  while(nrow(total_set[[length(total_set)]]) < s) {
    r <- runif(1)
    if (r < pm) {
      for (i in 1:M) {
        for (j in 1:length(total_set)) points <- mutate(points, all_specs, imp_func, imp_levels, index = j)
      }
    }
    else
      for (i in 1:floor((nrow(points)+1)/2)) points < crossover(points, imp_func, imp_levels)
    for (i in 1:(nrow(points)+1)) points <- exchange(points, imp_func, imp_levels)
    if (!all(points[nrow(points),] == start_points[nrow(points),])) {
      for (i in 1:length(total_set)) total_set[[i]] <- do.call('rbind', list(total_set[[i]], points[i,]))
      start_points <- points
    }
  }
  return(total_set)
}

#' IDEMC Point Generation
#'
#' Performs Implausibility-Driven Evolutionary Monte Carlo. Given a set of initial points
#' (preferably sampled across the full space), the implausibility ladder is set up via a burn-in
#' phase, before a full set of points is generated. This is a very compuationally intensive
#' procedure for generating points, and should be used only when the target space is expected
#' to be very small or have strange disconnected structure. For less awkward target spaces, use
#' any of the functionality in \code{\link{generate_new_runs}}.
#'
#' The burn-in starts with a rung defined as the full space (i.e. any point whose implausibility
#' is less than the maximum implausibility over the space); from the sample of points it finds the
#' value of the implausibility such that the proportion of points in the new rung is \code{p} times
#' the number in the previous rung. It then uses these points to generate \code{s} new points at
#' the new rung using \code{IDEMC}. This continues until a desired lower rung is found (defined
#' by points whose implausibility is lower than \code{imp}).
#'
#' Once the burn-in is performed, a full set of \code{sn} points is produced using this ladder.
#'
#' @param xsamp The initial points
#' @param ems The emulators to evaluate implausibility over
#' @param targets The corresponding output targets
#' @param s The number of points to generate at each burn-in stage
#' @param sn The final number of points to generate
#' @param p The proportion of points kept at each stage of burn-in
#' @param imp The value of implausibility to stop the ladder at
#' @param all_specs If burn-in has already been performed, the cluster specifications
#' @param imps If burn-in has alreadt been performed, the implausibility laddder values
#' @param ... Any additional parameters to pass to \code{IDEMC_step}
#'
#' @return A list of data.frames, corresponding to the points generated at each rung
#'
#' @seealso [generate_new_runs()] for other point generation mechanisms.
#'
#' @references
#' Vernon, I. & Williamson, D. (2013) Efficient uniform designs for multi-wave computer experiments. arXiv:1309.3520
#'
#' @examples
#' \dontrun{
#' ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#' out_vars <- c('nS', 'nI', 'nR')
#' o_ems <- emulator_from_data(GillespieSIR, out_vars, ranges)
#' t_ems <- purrr::map(seq_along(o_ems), ~o_ems[[.]]$adjust(GillespieSIR, out_vars[[.]]))
#' z <- list(
#'  nS = list(val = 281, sigma = 10.43),
#'  nI = list(val = 30, sigma = 11.16),
#'  nR = list(val = 689, sigma = 14.32)
#' )
#' start_pts <- data.frame(
#'   aSI = runif(500, ranges$aSI[1], ranges$aSI[2]),
#'   aIR = runif(500, ranges$aIR[1], ranges$aIR[2]),
#'   aSR = runif(500, ranges$aSR[1], ranges$aSR[2])
#' )
#' result <- IDEMC(start_pts, t_ems, z, 50, 100, 0.3, imp = 2)
#' }
#'
#' @export
#'
IDEMC <- function(xsamp, ems, targets, s, sn, p, imp = 3, all_specs = NULL, imps = NULL, ...) {
  test_i <- max(nth_implausible(ems, xsamp, targets, max_imp = Inf))
  check_imp <- function(x, imp) {
    for (i in 1:length(ems)) if (!ems[[i]]$implausibility(x, targets[[i]], imp)) return(FALSE)
    return(TRUE)
  }
  if (is.null(all_specs) || is.null(imps)) {
    imps <- c(test_i)
    specs_0 <- list(ranges = purrr::map(seq_along(xsamp[1,]), ~c(min(xsamp[,.]), max(xsamp[,.]))))
    sample_imps <- nth_implausible(ems, xsamp, targets, max_imp = Inf)
    samp_with_imps <- cbind(xsamp, sample_imps) %>% setNames(c(names(xsamp), "I"))
    samp_with_imps <- samp_with_imps[order(samp_with_imps$I),]
    test_index <- floor(p*nrow(samp_with_imps)) + 1
    test_i <- samp_with_imps[test_index, "I"]
    xsub <- samp_with_imps[1:test_index, names(xsamp)]
    all_specs <- list(specs_0, get_specs(xsub))
    chroms <- list(xsamp[nrow(xsamp),], xsub[nrow(xsub),])
    imps <- c(imps, test_i)
    while(test_i > imp) {
      new_sample <- IDEMC_step(do.call('rbind', chroms), check_imp, imps, all_specs, s, ...)
      new_sample <- new_sample[[length(new_sample)]]
      sample_imps <- nth_implausible(ems, new_sample, targets, max_imp = Inf)
      samp_with_imps <- cbind(new_sample, sample_imps) %>% setNames(c(names(new_sample), "I"))
      samp_with_imps <- samp_with_imps[order(samp_with_imps$I),]
      test_index <- floor(p*nrow(samp_with_imps))+1
      test_i <- samp_with_imps[test_index, "I"]
      if (test_i < imp) test_i <- imp
      xsub <- samp_with_imps[1:test_index, names(new_sample)]
      all_specs[[length(chroms)+1]] <- get_specs(xsub)
      chroms[[length(chroms)+1]] <- xsub[nrow(xsub),]
      imps <- c(imps, test_i)
    }
    print("Finished burn-in. Implausibility ladder:")
    print(imps)
  }
  else chroms <- xsamp
  final_out <- IDEMC_step(do.call('rbind', chroms), check_imp, imps, all_specs, sn, ...)
  for (i in 1:length(final_out)) row.names(final_out[[i]]) <- NULL
  return(final_out)
}

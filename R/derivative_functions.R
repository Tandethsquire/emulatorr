# A little helper function to get emulated derivative information at a point.
get_deriv_info <- function(em, x, ...) {
  deriv_exp <- purrr::map_dbl(names(em$ranges), function(y) em$get_exp(x, p = y, ...))
  deriv_var <- matrix(unlist(purrr::map(names(em$ranges), function(a) purrr::map_dbl(names(em$ranges), function(b) em$get_cov(x, p = a, pp = b, ...)))), nrow = length(em$ranges), byrow = T)
  return(list(exp = deriv_exp, var = deriv_var))
}


#' Derivative inner product (EXPERIMENTAL)
#'
#' Find the (uncertainty-modified) inner product between the derivative at a point \code{x} and
#' a proposal direction \code{v}.
#'
#' Given a point \code{x} and a direction \code{v}, we find the overlap between E[f'(x)]
#' and \code{v}. The emulated derivative has uncertainty
#' associated with it: this variance is given by v %*% Var[f'(x)] %*% v. Both of the
#' emulator quantities are obtained from \code{get_deriv_info}.
#'
#' This function is concerned with ascertaining whether a direction is oriented in the direction
#' of the emulator gradient. It allows for a consideration of 'emulated gradient descent'.
#'
#' @param em The emulator in question
#' @param x The point in input space from which we want to consider the derivative
#' @param v Any direction in the d-dimensional space considered by the emulator
#' @param ... Additional arguments to be passed along (e.g. local.var to the emulator)
#'
#' @return A 2-vector, consisting of the lower and upper (3-sigma) bounds for the inner product.
#'
#' @export
directional_fit <- function(em, x, v, ...) {
  normed_x <- x/sqrt(sum(x^2))
  normed_v <- v/sqrt(sum(v^2))
  deriv_info <- get_deriv_info(em, normed_x, ...)
  suitability <- deriv_info$exp %*% normed_v
  uncertainty <- sqrt(normed_v %*% deriv_info$var %*% normed_v)
  return(c(suitability - 3*uncertainty, suitability + 3*uncertainty))
}

#' Derivative Point Proposal (EXPERIMENTAL)
#'
#' Proposes new points based on old ones using derivative emulation.
#'
#' Given a point (preferably close to the implausibility boundary) \code{x}, we can
#' calculate the emulated gradient at the point for each emulator.
#' If the value of E[f(x)] is larger than the desired value, then this emulator wants
#' the point to travel in the negative gradient direction, and conversely for smaller
#' E[f(x)]. The combination of this information for each emulator defines a preferred
#' set of directions of travel from this point.
#'
#' We can try to find a shared direction which improves all emulator evaluations; if some
#' outputs are already well inside the implausibility cutoff (i.e. if their implausibility)
#' is less than \code{accept}, then we can allow these targets to get worse to make the others
#' better.
#'
#' Provided a shared direction v has been identified, we then move in this direction as follows.
#' The new point is defined as x + h*v, for some choice of h. To determine h, we initialise
#' h = 0.1, propose a new point, and check the new measure: the measures are either implausibility
#' or (default) the difference between the target value and the emulator expectation, normalised
#' by the target value. If the new measure
#' is lower than the original one, we step along the direction further. If it is not, we reduce
#' h and try again. This iteration terminates either when the new proposed point starts to
#' increase the implausibility, if the value of h is particularly small (as determined by
#' \code{hcutoff}), or if we have taken more than \code{iteration.steps} steps in the proposal.
#'
#' @param ems The list of \code{\link{Emulator}} objects
#' @param x The currently proposed point
#' @param targets The list of emulator targets
#' @param accept The implausibility value at which the proposal can increase the implausibility
#' @param hcutoff A power of 10: the smallest allowed step-size, h
#' @param iteration.measure Which measure to use for point suitability: expectation or implausibilty?
#' @param iteration.steps How many iterations to perform before returning the point.
#'
#' @return A new proposal point, or the original point if a suitable new point could not be found.
#'
#' @export
directional_proposal <- function(ems, x, targets, accept = 2, hcutoff = 1e-06, iteration.measure = 'exp', iteration.steps = 50) {
  point_implaus <- purrr::map_dbl(seq_along(ems), ~ems[[.]]$implausibility(x, z = targets[[.]]))
  x_predict <- purrr::map_dbl(ems, ~.$get_exp(x))
  is_bigger <- purrr::map_lgl(seq_along(targets), ~targets[[.]]$val < x_predict[[.]])
  x_diffs <- do.call('rbind', purrr::map(ems, function(y) purrr::map_dbl(names(x), ~y$get_exp(x, p = .))))
  x_dir <- do.call('rbind', purrr::map(seq_along(is_bigger), ~if(is_bigger[[.]]) -1*x_diffs[.,] else x_diffs[.,]))
  x_norms <- apply(x_dir, 1, function(y) sqrt(sum(y^2)))
  x_dir <- sweep(x_dir, 1, x_norms, "/")
  test_dirs <- runifs(500 * length(x), length(x))
  suits <- apply(x_dir, 1, function(y) apply(test_dirs, 1, function(z) z %*% y))
  suit_means <- apply(suits, 1, mean)
  order_dirs <- test_dirs[order(suit_means, decreasing = TRUE),]
  order_suits <- suits[order(suit_means, decreasing = TRUE),]
  restrict_dirs <- order_dirs[apply(order_suits, 1, function(y) all(y >= 0 || point_implaus < accept)),]
  nth_discrepancy <- function(ems, x, targets, n = 1) {
    discs <- purrr::map_dbl(seq_along(ems), ~abs((ems[[.]]$get_exp(x)-targets[[.]]$val)/targets[[.]]$val))
    if (n == 1) return(max(discs))
    else return(order(discs, decreasing = TRUE)[[n]])
  }
  if (nrow(restrict_dirs) == 0) {
    warning("No direction reduces all relevant targets.")
    return(x)
  }
  range_dists <- purrr::map_dbl(ems[[1]]$ranges, ~(.[[2]]-.[[1]])/2)
  best_dir <- restrict_dirs[1,] * range_dists
  old_measure <- if (iteration.measure == "exp") nth_discrepancy(ems, x, targets) else max(point_implaus)
  better_pt <- NULL
  attempts <- 0
  gap <- 0.1
  index <- 1
  while(attempts < 50) {
    new_point <- x + gap * index * best_dir
    new_measure <- if (iteration.measure == "exp") nth_discrepancy(ems, new_point, targets) else nth_implausible(ems, new_point, targets)
    if (new_measure < old_measure) {
      better_pt <- new_point
      old_measure <- new_measure
      index <- index + 1
    }
    else {
      if (is.null(better_pt) && gap > hcutoff) {
        gap <- gap * 0.1
      }
      else {
        break
      }
    }
    attempts <- attempts + 1
  }
  if (is.null(better_pt)) return(x)
  return(better_pt)
}

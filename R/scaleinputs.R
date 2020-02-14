#' Scale input parameter values
#'
#' Given ranges and a point x taken from an arbitrary range, scale it
#' to be in the range [-1,1]*n (n is the number of input dimensions).
#'
#' @param x The point to scale
#' @param r The list of ranges
#'
#' @return The scaled point
#' @export
#'
#' @examples
#'    x <- c(0.1, 0.5, 0.025)
#'    ranges <- list(c(0.1,0.8), c(0, 0.5), c(0, 0.05))
#'    scale_input(x, ranges)
scale_input <- function(x, r) {
  centers <- purrr::map_dbl(r, ~(.x[2]+.x[1])/2)
  scales <- purrr::map_dbl(r, ~(.x[2]-.x[1])/2)
  return((x-centers)/scales)
}

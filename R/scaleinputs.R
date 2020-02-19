#' Scale input parameter values
#'
#' Given ranges and a point x taken from an arbitrary range, scale it
#' to be in the range [-1,1]*n (n is the number of input dimensions).
#'
#' @param x The point to scale
#' @param r The list of ranges
#' @param forward scale to [-1,1] (or from [-1,1])?
#'
#' @return The scaled point
#' @export
#'
#' @examples
#'    x <- c(0.1, 0.5, 0.025)
#'    ranges <- list(c(0.1,0.8), c(0, 0.5), c(0, 0.05))
#'    scale_input(x, ranges)
#'
#'    df <- data.frame(x = c(0.1, 0.4, 0.56), y = c(0.2, 0.05, 0.45), c(0.01, 0.025, 0.04))
#'    df_scaled <- apply(df, 1, scale_input, ranges)
#'    apply(df_scaled, 2, scale_input, ranges, FALSE)
scale_input <- function(x, r, forward = TRUE) {
  centers <- purrr::map_dbl(r, ~(.x[2]+.x[1])/2)
  scales <- purrr::map_dbl(r, ~(.x[2]-.x[1])/2)
  if (forward)
   return((x-centers)/scales)
  return(x*scales + centers)
}

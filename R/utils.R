# Scales inputs: important since the emulators should take inputs purely in [-1,1]
scale_input <- function(x, r, forward = TRUE) {
  centers <- purrr::map_dbl(r, ~(.x[2]+.x[1])/2)
  scales <- purrr::map_dbl(r, ~(.x[2]-.x[1])/2)
  if (is.null(names(x))) {
    centers <- unlist(centers, use.names = F)
    scales <- unlist(scales, use.names = F)
  }
  if (forward)
    return((x-centers)/scales)
  return(x*scales + centers)
}

# Helper to convert functions to names
function_to_names <- function(f, var_names) {
  basis_vectors <- setNames(unique(expand.grid(split(diag(1, nrow = length(var_names)), rep(1:length(var_names), each=length(var_names))))),var_names)
  namegrid <- purrr::map(var_names, ~rep("",length(var_names)))
  for (i in 1:length(var_names)) namegrid[[i]][i] <- var_names[i]
  namegrid <- gsub("(^\\*|\\*$)", "", gsub("\\*+", "*", apply(unique(expand.grid(namegrid)), 1, paste, collapse="*")))
  ordering <- order(apply(basis_vectors,1,sum))
  basis_vectors <- basis_vectors[ordering,]
  namegrid <- c(namegrid[ordering], use.names = F)
  out_str <- ""
  for (i in 1:length(namegrid)) {
    if (f(basis_vectors[i,]) != 0)
    {
      return(ifelse(namegrid[i] == "", "1", unlist(namegrid[i])))
    }
  }
}

# Wrapper for evaluating multiple functions over multiple points
eval_funcs <- function(funcs, points, ...) {
  if (class(points) != "data.frame" && length(dim(points)) == 0) points <- matrix(points, nrow = length(points))
  if (typeof(funcs) == "closure")
    output <- apply(points, 1, funcs, ...)
  else
    output <- apply(points, 1, function(x) sapply(funcs, purrr::exec, x))
  if (length(dim(output)) != 0) return(t(output))
  return(output)
}

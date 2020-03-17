# Scales inputs: important since the emulators should take inputs purely in [-1,1]
scale_input <- function(x, r, forward = TRUE) {
  centers <- purrr::map(r, ~(.x[2]+.x[1])/2)
  scales <- purrr::map(r, ~(.x[2]-.x[1])/2)
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

# Evaluate multiple functions over points
eval_funcs <- function(funcs, points, ...) {
  output <- NULL
  if (typeof(funcs) == "closure")
  {
    try(output <- apply(points, 1, funcs, ...), silent = TRUE)
    if (is.null(output)) try(output <- purrr::exec(funcs, points, ...), silent = TRUE)
  }
  else
  {
    try(output <- apply(points, 1, function(x) purrr::map_dbl(funcs, purrr::exec, x, ...)), silent = TRUE)
    if(is.null(output)) try(output <- purrr::map_dbl(funcs, purrr::exec, points, ...), silent = TRUE)
  }
  if(is.null(output)) stop("Could not evaluate functions at points.")
  return(output)
}

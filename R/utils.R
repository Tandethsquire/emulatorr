# Implausibility colour palettes
redgreen <- c('#00FF00', '#18FF00', '#31FF00', '#49FF00', '#62FF00', '#7AFF00', '#93FF00', '#ABFF00', '#C4FF00', '#DDFF00',
              '#E0E200', '#E4C600', '#E8AA00', '#EC8D00', '#EF7100', '#F35500', '#F73800', '#FB1C00', '#FF0000', '#FF0000')
colourblind <- c('#1aff1a', '#2af219', '#3ae618', '#4ada17', '#5acd16', '#6bc115', '#7bb514', '#8ba813', '#9b9c12', '#ac9011',
                 '#a2831d', '#98762a', '#8e6936', '#845d43', '#7a504f', '#70435c', '#663768', '#5c2a75', '#521d81', '#48118e')
redgreencont <- list(low = '#00FF00', mid = '#DDFF00', high = '#FF0000')
colourblindcont <- list(low = '#1aff1a', mid = '#ac9011', high = '#48118e')

# Scales inputs: important since the emulators should take inputs purely in [-1,1]
scale_input <- function(x, r, forward = TRUE) {
  centers <- purrr::map(r, ~(.x[2]+.x[1])/2)
  scales <- purrr::map(r, ~(.x[2]-.x[1])/2)
  if (is.null(names(x))) {
    centers <- unlist(centers, use.names = F)
    scales <- unlist(scales, use.names = F)
  }
  if (forward)
    output <- (x-centers)/scales
  else
    output <- x*scales + centers
  if (!"data.frame" %in% class(output))
    return(data.frame(output))
  return(output)
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
  pointsdim <- (length(dim(points)) != 0)
  manyfuncs <- (typeof(funcs) != "closure")
  if (manyfuncs && pointsdim)
    return(apply(points, 1, function(x) purrr::map_dbl(funcs, purrr::exec, x, ...)))
  if (manyfuncs)
    return(purrr::map_dbl(funcs, purrr::exec, points, ...))
  if (pointsdim) {
    return(tryCatch(apply(points, 1, funcs, ...),
                    error = function(cond1) {
                      tryCatch(purrr::exec(funcs, points, ...),
                               error = function(cond2) {
                                 print(cond1, cond2)
                               })
                    }))
  }
  return(purrr::exec(funcs, points, ...))
}

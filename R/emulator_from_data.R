#' Create Bayes Linear emulators from data
#'
#' Given simulator data, generates a set of emulators.
#'
#' Many of the parameters passed to the function are optional; if they are not
#' provided, a linear model is fitted as follows:
#'
#' The basis functions are assumed to be (1, x1, x2, ...) for inputs x1, x2, ...
#'
#' The associated beta parameters for each output are found from the coefficients of
#' the model fits from \code{\link{get_linear_model}} or \code{\link{get_quadratic_model}},
#' depending on the value of the boolean \code{quadratic}.
#'
#' The correlation structure is taken to be exponential squared, with correlation length one-sixth
#' of the input parameter range. The overall variance is taken to be the residual squared error
#' from the linear models used to fit the regression coefficients, and zero expectation is assumed.
#' If a nugget \code{w(x)} is desired, then the scale of the nugget should be provided as a list
#' of numerics to \code{deltas} (see the \code{\link{Correlator}} documentation for details).
#'
#' The covariance between the regression coefficients and the correlation structure is assumed to
#' be 0 (passed as \code{NULL} into the emulator).
#'
#' @param input_data Required. A dataframe containing the input parameters and output values from a set of simulator runs.
#' @param input_names Required. Indicates which columns of \code{input_data} are the inputs.
#' @param output_names Required. Indicates which columns of \code{input_data} are the outputs.
#' @param ranges The input parameter ranges, as a list of 2-vectors.
#' @param beta Optional. Specifications for regression coefficients; each in the form \code{list(mu = vector, sigma = matrix)} per output.
#' @param u Optional. Specifications for correlation structure; each in the form \code{list(sigma = numeric, mu = closure, theta = numeric, corr = closure)}
#' @param c_lengths Optional if \code{ranges} or a full \code{u} specification is provided. Correlation lengths for each output.
#' @param funcs Optional. Basis functions for the regression surface.
#' @param bucov Optional. Covariance vector between beta and u.
#' @param quadratic Include quadratic terms in model fitting? Default: FALSE
#' @param deltas Optional. Specifications for the nugget terms.
#'
#' @return A list of objects of class \code{\link{Emulator}}.
#' @export
#'
#' @examples
#'     inputdata <- GillespieSIR
#'     out_vars <- c("nS", "nI", "nR")
#'     ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#'     emulators <- emulator_from_data(input_data = inputdata, input_names = names(ranges),
#'      output_names = out_vars, ranges = ranges, c_lengths = c(3/4,3/4,3/4))
#'     emulators[[1]]$get_exp(c(0.4,0.25,0.025))
#'     #> 640.9275
#'     emulators[[1]]$get_cov(c(0.4,0.25,0.025))
#'     #> 8939.337
#'     # Now we can actually use the data to generate useful emulators
#'     new_emulators <- purrr::map(seq_along(emulators),
#'         ~emulators[[.x]]$bayes_adjust(inputdata[,names(ranges)], inputdata[,out_vars[[.x]]]))
#'     new_emulators[[1]]$get_exp(c(0.4,0.25,0.025))
#'     #> 344.4144
#'     new_emulators[[1]]$get_cov(c(0.4,0.25,0.025))
#'     #> 8902.77
#'
#'     # Same data as above, but with custom specifications.
#'     # Suppose we instead use GLS estimates of the betas, and they're known:
#'     beta_sigma <- diag(0, nrow=4)
#'     beta_specs <- list(
#'      list(mu = c(508.8783, -616.2146, 863.3770, -5616.7661), sigma = beta_sigma),
#'      list(mu = c(287.0726, 251.0185, -1029.0808, -705.5883), sigma = beta_sigma),
#'      list(mu = c(183.7949, 429.8564, 191.6515, 6058.9556), sigma = beta_sigma)
#'     )
#'     # Generate the correlators
#'     u_mu <- function(x) 0
#'     correlators <- list(
#'      list(sigma = 94.548, mu = u_mu, theta = 3/4, corr = exp_sq),
#'      list(sigma = 113.299, mu = u_mu, theta = 3/4, corr = exp_sq),
#'      list(sigma = 144.198, mu = u_mu, theta = 3/4, corr = exp_sq)
#'     )
#'     basis_functions <- c(function(x) 1, function(x) x[[1]],
#'      function(x) x[[2]], function(x) x[[3]])
#'     # Use nugget terms based on sd of training data
#'     deltas <- c(0.08522, 0.00604, 0.04123)
#'     emulators <- emulator_from_data(input_data = GillespieSIR,
#'      input_names = names(ranges),
#'      output_names = out_vars, beta = beta_specs, u = correlators,
#'      funcs = basis_functions, ranges = ranges, deltas = deltas)
#'     emulators[[1]]$get_exp(c(0.4, 0.25, 0.025))
#'     #> 596.909
#'     emulators[[1]]$get_cov(c(0.4, 0.25, 0.025))
#'     #> 8939.324
#'
#'     # Alternatively, allow quadratic pieces:
#'     quadratic_emulators <- emulator_from_data(input_data = inputdata,
#'      input_names = names(ranges), ranges = ranges,
#'      output_names = out_vars, quadratic = TRUE, deltas = deltas)
#'     quadratic_emulators[[1]]$get_exp(c(0.4, 0.25, 0.025))
#'     #> 339.8135
#'     quadratic_emulators[[1]]$get_cov(c(0.4, 0.25, 0.025))
#'     #>  3890.263
#'
emulator_from_data <- function(input_data, input_names, output_names, ranges, beta, u, c_lengths, funcs, bucov, deltas, quadratic=F) {
  if (missing(ranges)) {
    warning("No ranges provided; inputs assumed to be in the range [-1,1].")
    ranges <- purrr::map(input_names, ~c(-1,1))
  }
  data <- cbind(t(apply(input_data[,input_names], 1, scale_input, ranges)), input_data[,output_names])
  if (missing(beta) || missing(u) || missing(funcs)) {
    if (quadratic) {
      does_add <- choose(length(input_names)+2, length(input_names))>length(input_data[,1])
      models <- purrr::map(output_names, ~get_quadratic_model(input_data, ranges, .x, add = does_add))
    }
    else
    {
      does_add <- length(input_names)+1>length(input_data[,1])
      models <- purrr::map(output_names, ~get_linear_model(input_data, ranges, .x, add = does_add))
    }
  }
  if (missing(beta)) {
    all_funcs <- c(function(x) 1, purrr::map(1:length(input_names), ~function(x) x[[.x]]))
    all_coeffs <- c("(Intercept)", input_names)
    if (quadratic) {
      all_funcs <- c(all_funcs, apply(expand.grid(1:length(input_names), 1:length(input_names)), 1, function(y) function(x) x[[y[[1]]]]*x[[y[[2]]]]))
      all_coeffs <- c(all_coeffs, apply(expand.grid(input_names, input_names), 1, paste, collapse=":"))
    }
    model_beta_mus <- purrr::map(models, ~c(.x$coefficients, use.names = F))
    model_basis_funcs <- purrr::map(models, ~all_funcs[all_coeffs %in% variable.names(.x)])
    model_beta_sigmas <- purrr::map(model_beta_mus, ~diag(0, nrow = length(.x)))
  }
  else {
    model_beta_mus <- purrr::map(beta, ~.x$mu)
    model_beta_sigmas <- purrr::map(beta, ~.x$sigma)
    model_basis_funcs <- purrr::map(beta, ~funcs)
  }
  if (missing(u)) {
    model_u_sigmas <- c(sapply(models, function(x) summary(x)$sigma))
    model_u_mus <- purrr::map(output_names, ~function(x) 0)
    if(missing(c_lengths)) {
      ifelse(quadratic, model_u_thetas <- runif(length(output_names), 1/3, 2/3), model_u_thetas <- runif(length(output_names), 1/2, 1))
    }
    else model_u_thetas <- c_lengths
    model_u_corr_funcs <- purrr::map(output_names, ~exp_sq)
  }
  else {
    model_u_sigmas <- purrr::map(u, ~.x$sigma)
    model_u_mus <- purrr::map(u, ~.x$mu)
    model_u_thetas <- purrr::map(u, ~.x$theta)
    model_u_corr_funcs <- purrr::map(u, ~.x$corr)
  }
  if (missing(deltas))
    model_us <- purrr::map(seq_along(model_u_sigmas), ~Correlator$new(function(x, y) model_u_sigmas[[.x]]^2*model_u_corr_funcs[[.x]](x, y, model_u_thetas[[.x]]), model_u_mus[[.x]]))
  else
    model_us <- purrr::map(seq_along(model_u_sigmas), ~Correlator$new(function(x, y) model_u_sigmas[[.x]]^2*model_u_corr_funcs[[.x]](x, y, model_u_thetas[[.x]]), model_u_mus[[.x]], deltas[[.x]]))
  model_betas <- purrr::map(seq_along(model_beta_mus), ~list(mu = model_beta_mus[[.x]], sigma = model_beta_sigmas[[.x]]))
  if (missing(bucov))
    out_ems <- purrr::map(seq_along(model_betas), ~Emulator$new(funcs = model_basis_funcs[[.x]], beta = model_betas[[.x]], u = model_us[[.x]], ranges = ranges))
  else
    out_ems <- purrr::map(seq_along(model_betas), ~Emulator$new(funcs = model_basis_funcs[[.x]], beta = model_betas[[.x]], u = model_us[[.x]], bucov = bucov, ranges = ranges))
  return(out_ems)
}

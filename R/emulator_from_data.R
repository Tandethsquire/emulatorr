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
#' step(lm(output~x1+x2+...)).
#' They are assumed to be known, and so given a variance matrix 0.
#'
#' The correlation structire is taken to be exponential squared, with correlation length one-sixth
#' of the input parameter range. The overall variance is taken to be the residual squared error
#' from the linear models used to fit the regression coefficients, and zero expectation is assumed.
#'
#' The covariance between the regression coefficients and the correlation structure is assumed to
#' be 0 (passed as \code{NULL} into the emulator).
#'
#' @importFrom purrr %>%
#' @importFrom stats as.formula lm step
#'
#' @param input_data Required. A dataframe containing the input parameters and output values from a set of simulator runs.
#' @param input_names Required. Indicates which columns of \code{input_data} are the inputs.
#' @param output_names Required. Indicates which columns of \code{input_data} are the outputs.
#' @param ranges Optional if \code{c_lengths} is provided. The input parameter ranges, as a list of 2-vectors.
#' @param beta Optional. Specifications for regression coefficients; each in the form \code{list(mu = vector, sigma = matrix)} per output.
#' @param u Optional. Specifications for correlation structure; each in the form \code{list(sigma = numeric, mu = closure, theta = numeric, corr = closure)}
#' @param c_lengths Optional if \code{ranges} or a full \code{u} specification is provided. Correlation lengths for each output.
#' @param funcs Optional. Basis functions for the regression surface.
#' @param bucov Optional. Covariance vector between beta and u.
#'
#' @return A list of objects of class \code{\link{Emulator}}.
#' @export
#'
#' @examples
#'     inputdata <- GillespieSIR
#'     in_vars <- c("aSI", "aIR", "aSR")
#'     out_vars <- c("nS", "nI", "nR")
#'     ranges <- list(c(0.1, 0.8), c(0, 0.5), c(0, 0.05))
#'     emulators <- emulator_from_data(input_data = inputdata, input_names = in_vars,
#'      output_names = out_vars, ranges = ranges, c_lengths = c(0.1, 0.085, 0.075))
#'     emulators[[1]]$get_exp(c(0.4,0.25,0.025))
#'     #> 341.8855
#'     emulators[[1]]$get_var(c(0.4,0.25,0.025))
#'     #> 8939.337
#'     # Now we can actually use the data to generate useful emulators
#'     new_emulators <- purrr::map(seq_along(emulators),
#'         ~emulators[[.x]]$bayes_adjust(inputdata[,in_vars], inputdata[,out_vars[[.x]]]))
#'     new_emulators[[1]]$get_exp(c(0.4,0.25,0.025))
#'     #> 367.8078
#'     new_emulators[[1]]$get_var(c(0.4,0.25,0.025))
#'     #> 679.3054
#'
#'     # Same data as above, but with custom specifications.
#'     # Suppose we instead use GLS estimates of the betas, and they're known:
#'     beta_sigma <- diag(0, nrow=4)
#'     beta_specs <- list(
#'      list(mu = c(508.8783, -616.2146, 863.3770, -5616.7661), sigma = beta_sigma),
#'      list(mu = c(287.0726, 251.0185, -1029.0808, -705.5883), sigma = beta_sigma),
#'      list(mu = c(183.7949, 429.8564, 191.6515, 6058.9556), sigma = beta_sigma)
#'     )
#'     u_mu <- function(x) 0
#'     correlators <- list(
#'      list(sigma = 94.548, mu = u_mu, theta = 0.1, corr = exp_sq),
#'      list(sigma = 113.299, mu = u_mu, theta = 0.085, corr = exp_sq),
#'      list(sigma = 144.198, mu = u_mu, theta = 0.075, corr = exp_sq)
#'     )
#'     emulators <- emulator_from_data(input_data = GillespieSIR, input_names = in_vars,
#'      output_names = out_vars, beta = beta_specs, u = correlators)
#'
emulator_from_data <- function(input_data, input_names, output_names, ranges, beta, u, c_lengths, funcs, bucov) {
  if (missing(beta) || missing(u) || missing(funcs)) {
    model_params <- paste(input_names, collapse = "+")
    models <- input_data[output_names] %>% apply(MARGIN = 2, FUN = function(x) step(lm(as.formula(paste('x~', model_params, sep = "")), data = input_data), trace = 0))
  }
  if (missing(beta)) {
    coeffnames <- c('(Intercept)', input_names)
    model_beta_mu <- purrr::map(models, ~ifelse(is.na(.x$coefficients[coeffnames]), 0, .x$coefficients[coeffnames] %>% as.numeric))
    model_beta_mu <- purrr::map(seq_along(model_beta_mu), ~c(model_beta_mu[.x][[1]] %>% as.numeric))
    model_beta_sigma <- purrr::map(seq_along(1:length(input_names)), ~diag(0, nrow=length(model_beta_mu[[1]])))
  }
  else {
    model_beta_mu <- purrr::map(seq_along(beta), ~beta[[.x]]$mu)
    model_beta_sigma <- purrr::map(seq_along(beta), ~beta[[.x]]$sigma)
  }
  if (missing(u)) {
    model_u_sigma <- c(sapply(models, function(x) summary(x)$sigma))
    model_u_mu <- purrr::map(seq_along(1:length(input_names)), .f=~function(x) 0)
    if (missing(c_lengths)) model_u_theta <- purrr::map_dbl(ranges, ~(.x[2]-.x[1])/6)
    else model_u_theta <- c_lengths
    model_u_corr_func <- purrr::map(seq_along(1:length(input_names)), ~exp_sq)
  }
  else {
    model_u_sigma <- purrr::map(seq_along(u), ~u[[.x]]$sigma)
    model_u_mu <- purrr::map(seq_along(u), ~u[[.x]]$mu)
    model_u_theta <- purrr::map(seq_along(u), ~u[[.x]]$theta)
    model_u_corr_func <- purrr::map(seq_along(u), ~u[[.x]]$corr)
  }
  if (missing(funcs)) {
    model_funcs <- c(function(x) 1, purrr::map(seq_along(1:length(input_names)), ~function(x) x[[.x]]))
  }
  else {
    model_funcs <- funcs
  }
  model_u <- purrr::map(seq_along(model_beta_mu), ~list(theta = model_u_theta[[.x]], sigma = model_u_sigma[[.x]], c_func = model_u_corr_func[[.x]], mu_func = model_u_mu[[.x]]))
  model_beta <- purrr::map(seq_along(model_beta_mu), ~list(mu = model_beta_mu[[.x]], sigma = model_beta_sigma[[.x]]))
  correlators <- purrr::map(model_u, ~Correlator$new(function(x, y) .x$sigma^2*.x$c_func(x, y, .x$theta), .x$mu_func))
  if (missing(bucov))
    output_emulators <- purrr::map(seq_along(model_beta), ~Emulator$new(model_funcs, list(mu = model_beta[[.x]]$mu, sigma = model_beta[[.x]]$sigma), correlators[[.x]]))
  else
    output_emulators <- purrr::map(seq_along(model_beta), ~Emulator$new(model_funcs, list(mu = model_beta[[.x]]$mu, sigma = model_beta[[.x]]$sigma), correlators[[.x]], bucov))
  return(output_emulators)
}

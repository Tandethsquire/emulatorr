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
#' @param ranges The input parameter ranges, as a list of 2-vectors.
#' @param beta Optional. Specifications for regression coefficients; each in the form \code{list(mu = vector, sigma = matrix)} per output.
#' @param u Optional. Specifications for correlation structure; each in the form \code{list(sigma = numeric, mu = closure, theta = numeric, corr = closure)}
#' @param c_lengths Optional if \code{ranges} or a full \code{u} specification is provided. Correlation lengths for each output.
#' @param funcs Optional. Basis functions for the regression surface.
#' @param bucov Optional. Covariance vector between beta and u.
#' @param quadratic Include quadratic terms in model fitting? Default: FALSE
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
#'     basis_functions <- c(function(x) 1, function(x) x[[1]],
#'      function(x) x[[2]], function(x) x[[3]])
#'     emulators <- emulator_from_data(input_data = GillespieSIR, input_names = in_vars,
#'      output_names = out_vars, beta = beta_specs, u = correlators,
#'      funcs = basis_functions, ranges = ranges)
#'
emulator_from_data <- function(input_data, input_names, output_names, ranges, beta, u, c_lengths, funcs, bucov, quadratic=F) {
  if (missing(ranges)) ranges <- purrr::map(input_names, ~c(-1,1))
  data <- cbind(t(apply(input_data[,input_names], 1, scale_input, ranges)), input_data[,output_names])
  if (missing(beta) || missing(u) || missing(funcs)) {
    if(quadratic) form <- paste('(', paste(input_names, collapse="+"), ")^2", sep="")
    else form <- paste(input_names, collapse="+")
    models <- input_data[,output_names] %>% apply(MARGIN = 2, FUN = function(x) step(lm(as.formula(paste('x~',form,sep="")), data=input_data), trace = 0))
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
    if(missing(c_lengths)) model_u_thetas <- purrr::map_dbl(output_names, ~1/3)
    else model_u_thetas <- c_lengths
    model_u_corr_funcs <- purrr::map(output_names, ~exp_sq)
  }
  else {
    model_u_sigmas <- purrr::map(u, ~.x$sigma)
    model_u_mus <- purrr::map(u, ~.x$mu)
    model_u_thetas <- purrr::map(u, ~.x$theta)
    model_u_corr_funcs <- purrr::map(u, ~.x$corr)
  }
  model_us <- purrr::map(seq_along(model_u_sigmas), ~Correlator$new(function(x, y) model_u_sigmas[[.x]]^2*model_u_corr_funcs[[.x]](x, y, model_u_thetas[[.x]]), model_u_mus[[.x]]))
  model_betas <- purrr::map(seq_along(model_beta_mus), ~list(mu = model_beta_mus[[.x]], sigma = model_beta_sigmas[[.x]]))
  if (missing(bucov))
    out_ems <- purrr::map(seq_along(model_betas), ~Emulator$new(funcs = model_basis_funcs[[.x]], beta = model_betas[[.x]], u = model_us[[.x]], ranges = ranges))
  else
    out_ems <- purrr::map(seq_along(model_betas), ~Emulator$new(funcs = model_basis_funcs[[.x]], beta = model_betas[[.x]], u = model_us[[.x]], bucov = bucov, ranges = ranges))
  return(out_ems)
}

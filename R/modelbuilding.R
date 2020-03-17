#' Linear Model Generation
#'
#' Creates a linear model from data and a list of inputs.
#' There are two ways to generate the model; either start with all possible linear
#' terms, and then stepwise remove them (using \code{step}), or start with an intercept
#' add stepwise add linear terms, only retaining them if the AIC is improved. Which is
#' chosen is dependent on the value of \code{add}; in the event where \code{add = FALSE}
#' and there are not enough degrees of freedom, a warning will be given.
#'
#' @importFrom stats lm step setNames AIC update as.formula
#'
#' @param data A \code{data.frame} containing the input and output values
#' @param ranges A named \code{list} consisting of the ranges of the input parameters
#' @param output_name A string corresponding to the output to be modelled
#' @param add Should we perform stepwise add or stepwise delete? Default: \code{FALSE}
#'
#' @return The fitted model
#' @export
get_linear_model <- function(data, ranges, output_name, add = FALSE) {
  if (!add & (length(ranges) + 1 > length(data[,1]))) stop("Number of linear terms is greater than the degrees of freedom in the data. Consider using add = TRUE.")
  input_params <- names(ranges)
  scaled_input_data <- eval_funcs(scale_input, data[,names(ranges)], ranges)
  full_scaled <- setNames(data.frame(cbind(scaled_input_data, data[,output_name])), c(names(ranges), output_name))
  if (!add) {
    model <- step(lm(data = full_scaled, as.formula(paste(output_name, paste(names(ranges), collapse = "+"), sep="~"))), trace = 0)
  }
  else {
    model <- lm(data = full_scaled, as.formula(paste(output_name,  " ~ 1", sep = "")))
    aic <- AIC(model)
    for (i in 1:length(ranges)) {
      temp_model <- update(model, as.formula(paste(". ~ . + ", names(ranges)[i], sep = "")))
      temp_aic <- AIC(temp_model)
      if (temp_aic < aic) {
        model <- temp_model
        aic <- temp_aic
      }
    }
  }
  return(model)
}

#' Quadratic Model Generation
#'
#' Creates a quadratic model from data and a list of inputs.
#' There are two ways to generate the model; either start with all possible linear and
#' quadratic terms, and then stepwise remove them (using \code{step}), or start with a
#' linear model (maybe from \code{\link{get_linear_model}}) and add quadratic terms one
#' by one, only retaining them if the AIC is improved. Which is chosen is dependent on
#' the value of \code{add}; in the event where \code{add = FALSE} and there are not
#' enough degrees of freedom, a warning will be given.
#'
#' @importFrom stats lm step setNames AIC update as.formula
#'
#' @param data A \code{data.frame} containing the input and output values
#' @param ranges A named \code{list} consisting of the ranges of the input parameters
#' @param output_name A string corresponding to the output to be modelled
#' @param add Should we perform stepwise add or stepwise delete? Default: \code{FALSE}
#' @param linear_model Optional. A linear model to augment if add = TRUE
#'
#' @return The fitted model
#' @export
get_quadratic_model <- function(data, ranges, output_name, add = FALSE, linear_model = NULL) {
  if (!add & (choose(length(ranges)+2, length(ranges)) > length(data[,1]))) stop("Number of regression terms is greater than the degrees of freedom in the data. Use add = TRUE.")
  input_params <- names(ranges)
  scaled_input_data <- eval_funcs(scale_input, data[,names(ranges)], ranges)
  full_scaled <- setNames(data.frame(cbind(scaled_input_data, data[,output_name])), c(names(ranges), output_name))
  if (!add) {
    model <- step(lm(data = full_scaled, as.formula(paste(c(output_name, "~(", paste(names(ranges), collapse = "+"), ")^2"), collapse=""))), trace = 0)
  }
  else {
    ifelse (is.null(linear_model), model <- get_linear_model(data, ranges, output_name, add = add), model <- linear_model)
    aic <- AIC(model)
    poss_q <- apply(expand.grid(names(ranges), names(ranges)), 1, paste, collapse = ":")
    for (i in 1:length(poss_q)) {
      temp_model <- update(model, as.formula(paste(". ~ . + ", poss_q[i], sep = "")))
      temp_aic <- AIC(temp_model)
      if (temp_aic < aic) {
        model <- temp_model
        aic <- temp_aic
      }
    }
  }
  return(model)
}

#' Generate Prior Emulators from Data
#'
#' Given data from a simulation, generates a set of \code{\link{Emulator}} objects based on
#' fitted values.
#'
#' Many of the parameters that can be passed to this function are optional; the bare minimum
#' is \code{input_data}, \code{output_names}, and one of \code{ranges} or \code{input_names}.
#' If \code{ranges} is specified, then the input names are taken from that; if only
#' \code{input_names} is specified, then it is assumed that all input values in \code{input_data}
#' are already scaled to [-1, 1].
#'
#' If the minimum information is provided, then a model is fitted as follows.
#'
#' The basis functons and regression coefficients are generated using the \code{lm} function using
#' either only linear terms or up to quadratic terms (dependent on the value of \code{quadratic}),
#' performing stepwise add or delete as appropriate; in either event, the AIC criteria is used to
#' select the terms. The regression parameters thus derived are assumed to be known, so that
#' \code{beta$sigma = diag(0)}.
#'
#' The correlation function c(x,x') is taken to be \code{\link{exp_sq}}; the correlation length is
#' chosen using the Durham heuristic: this states that the correlation length should lie within
#' [1/(n+1), 2/(n+1)] where n is the degree of the fitted surface (and the range of the parameter
#' is [-1,1]). A value in this range is sampled uniformly: it may be useful to generate emulators
#' allowing this randomness and use diagnostics to decide on a fixed value. The nugget is assumed
#' to be 0 unless otherwise specified. The expectation E[u(x)] is assumed to be 0, and the
#' variance \code{sigma^2} is taken from the residual squared error from the model used to fit the
#' basis functions and betas.
#'
#' The covariance between beta and u(x) is assumed to vanish.
#'
#' @param input_data Required. A \code{data.frame} containing the input parameters and output values
#' from a set of simulator runs.
#' @param input_names A list of input_names (if \code{ranges} is not provided).
#' @param output_names Required. The list of outputs to emulate from \code{input_data}.
#' @param ranges A named list of parameter ranges.
#' @param beta Optional: specifications for the regression coefficients, given as a list
#' of lists \code{list(mu, sigma)} (a la \code{\link{Emulator}} specification).
#' @param u Optional: the correlation structure for each output, given as a list of
#' lists \code{list(mu, sigma, corr)}.
#' @param c_lengths Optional: a set of correlation lengths.
#' @param funcs Optional: basis functions for the regression surface.
#' @param bucov Optional: a list of functions giving the covariance between each of the beta
#' parameters and u(x).
#' @param deltas Optional: the nugget terms to include in u(x).
#' @param quadratic Optional: should the regression surface be linear or quadratic? Default: F
#'
#' @return A list of objects of class \code{\link{Emulator}}.
#' @export
#'
#' @examples
#'  # Use the GillespieSIR dataset
#'  ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#'  out_vars <- c('nS', 'nI', 'nR')
#'  ems_linear <- emulator_from_data(GillespieSIR, output_names = out_vars,
#'   ranges = ranges)
#'  ems_linear # Printout of the key information
#'
#'  ems <- emulator_from_data(GillespieSIR, output_names = out_vars,
#'   ranges = ranges, quadratic = TRUE)
#'  ems # Now includes quadratic terms (but only where they're warranted)
#'
#'  ems2 <- emulator_from_data(GillespieSIR, output_names = out_vars,
#'   ranges = ranges, c_lengths = c(0.55, 0.6, 0.59),
#'   deltas = c(0.1, 0.2, 0.2), quadratic = TRUE)
#'  ems2 # Broadly the same, but with the correlation structure modified.
emulator_from_data <- function(input_data, input_names = names(ranges), output_names, ranges, beta, u, c_lengths, funcs, bucov, deltas, quadratic=F) {
  if (missing(ranges)) {
    warning("No ranges provided; inputs assumed to be in the range [-1,1].")
    ranges <- purrr::map(input_names, ~c(-1,1))
  }
  data <- cbind(eval_funcs(scale_input, input_data[,names(ranges)], ranges), input_data[,output_names])
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
    model_u_corr_funcs <- purrr::map(seq_along(output_names), ~function(x, xp) exp_sq(x, xp, model_u_thetas[[.x]]))
  }
  else {
    model_u_sigmas <- purrr::map(u, ~.x$sigma)
    model_u_mus <- purrr::map(u, ~.x$mu)
    model_u_thetas <- purrr::map(u, ~.x$theta)
    model_u_corr_funcs <- purrr::map(u, ~.x$corr)
  }
  model_us <- purrr::map(seq_along(model_u_sigmas), ~list(mu = model_u_mus[[.x]], sigma = model_u_sigmas[[.x]], corr = model_u_corr_funcs[[.x]]))
  model_betas <- purrr::map(seq_along(model_beta_mus), ~list(mu = model_beta_mus[[.x]], sigma = model_beta_sigmas[[.x]]))
  if (missing(deltas)) model_deltas <- rep(0, length(output_names))
  else model_deltas <- deltas
  if (missing(bucov))
    out_ems <- purrr::map(seq_along(model_betas), ~Emulator$new(basis_f = model_basis_funcs[[.x]], beta = model_betas[[.x]], u = model_us[[.x]], ranges = ranges, delta = model_deltas[[.x]]))
  else
    out_ems <- purrr::map(seq_along(model_betas), ~Emulator$new(basis_f = model_basis_funcs[[.x]], beta = model_betas[[.x]], u = model_us[[.x]], bucov = bucov, ranges = ranges, delta = model_deltas[[.x]]))
  return(out_ems)
}

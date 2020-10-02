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
  if (!add & (length(ranges) + 1 > length(data[,1]))) {
    warning("Number of linear terms is greater than the degrees of freedom in the data. Changing to add = TRUE.")
    add <- TRUE
  }
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
#' @import purrr
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
  .x <- NULL
  if (!add & (choose(length(ranges)+2, length(ranges)) > length(data[,1]))) {
    warning("Number of regression terms is greater than the degrees of freedom in the data. Changing to add = TRUE.")
    add <- TRUE
  }
  input_params <- names(ranges)
  scaled_input_data <- eval_funcs(scale_input, data[,names(ranges)], ranges)
  full_scaled <- setNames(data.frame(cbind(scaled_input_data, data[,output_name])), c(names(ranges), output_name))
  if (!add) {
    quad_terms <- paste(purrr::map_chr(names(ranges), ~paste("I(",.x,"^2)", sep = "")), collapse = "+", sep = "")
    model <- step(lm(data = full_scaled, as.formula(paste(c(output_name, "~(", paste(names(ranges), collapse = "+"), ")^2+", quad_terms), collapse=""))), trace = 0)
  }
  else {
    ifelse (is.null(linear_model), model <- get_linear_model(data, ranges, output_name, add = add), model <- linear_model)
    aic <- AIC(model)
    quad_terms <- purrr::map_chr(names(ranges), paste("I(",.x,"^2)", sep =  ""))
    poss_q <- c(apply(expand.grid(names(ranges), names(ranges)), 1, paste, collapse = ":"), quad_terms)
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

# Performs maximum likelihood estimation of the hyperparameters sigma and theta.
# At present, only does single theta estimates (i.e. the same correlation length
# is applied in all parameter directions)
get_likelihood <- function(inputs, outputs, h, theta_range, n_ints = 200, delta = 0.1) {
  best_L <- -Inf
  best_sigma <- best_theta <- NULL
  hs <- eval_funcs(h, inputs)
  corr_func <- function(x, xp, theta, delta) {
    (1-delta) * exp(-sum((x-xp)^2/theta^2)) + ifelse(all(x == xp), delta, 0)
  }
  corr_matrix <- function(points, theta, delta) {
    apply(points, 1, function(x) apply(points, 1, function(y) corr_func(x,y,theta,delta)))
  }
  t_list <- seq(theta_range[1], theta_range[2], length.out = n_ints)
  for (i in t_list) {
    cmat <- corr_matrix(inputs, i, delta)
    cm <- chol2inv(chol(cmat))
    #beta_est <- chol2inv(chol(hs %*% cm %*% t(hs))) %*% hs %*% cm %*% outputs
    beta_est <- rep(0, nrow(hs))
    var_est <- t(outputs - t(hs) %*% beta_est) %*% cm %*% (outputs - t(hs) %*% beta_est)
    lik <- det(cmat)^(-1/2) * det(hs %*% cm %*% t(hs))^(-1/2) * var_est^(-(nrow(inputs)-nrow(hs))/2)
    if (lik > best_L) {
      best_L <- lik
      best_sigma <- sqrt(var_est)
      best_theta <- i
    }
  }
  return(list(sigma = best_sigma/sqrt(length(outputs)), theta = best_theta))
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
#' select the terms. The regression parameters thus derived are assumed to be known if
#' \code{beta.var=FALSE}, so that \code{beta$sigma = diag(0)}. Otherwise, the covariance matrix
#' for the parameters is taken from \code{vcov(model)}.
#'
#' The correlation function c(x,x') is taken to be \code{\link{exp_sq}}; the correlation length is
#' chosen using the Durham heuristic: this states that the correlation length should lie within
#' [1/(n+1), 2/(n+1)] where n is the degree of the fitted surface (and the range of the parameter
#' is [-1,1]). Maximum likelihod estimation is then applied to this range to find an acceptable
#' correlation length, and the corresponding standard error is used as an estimate for the variance
#' of the correlation structure. The expectation E[u(x)] is assumed to be 0.
#'
#' If delta terms are provided, then the nugget terms for each emulator are defined using these.
#' If they are not provided but a list of variabilities for each output are (in \code{ev}), then
#' a rough estimate of the nugget terms is performed and the emulators obtain these terms. If
#' neither is provided, the nugget terms are assumed to be identically zero for each emulator.
#'
#' The covariance between beta and u(x) is assumed to vanish.
#'
#' @param input_data Required. A \code{data.frame} containing the input parameters and output values
#' from a set of simulator runs.
#' @param output_names Required. The list of outputs to emulate from \code{input_data}.
#' @param ranges A named list of parameter ranges.
#' @param input_names A list of input_names (if \code{ranges} is not provided).
#' @param beta Optional: specifications for the regression coefficients, given as a list
#' of lists \code{list(mu, sigma)} (a la \code{\link{Emulator}} specification).
#' @param u Optional: the correlation structure for each output, given as a list of
#' lists \code{list(mu, sigma, corr)}.
#' @param c_lengths Optional: a set of correlation lengths.
#' @param funcs Optional: basis functions for the regression surface.
#' @param bucov Optional: a list of functions giving the covariance between each of the beta
#' parameters and u(x).
#' @param deltas Optional: the nugget terms to include in u(x).
#' @param ev Optional. Used for determining nugget terms in absence on delta
#' @param quadratic Optional: should the regression surface be linear or quadratic? Default: F
#' @param beta.var Optional: should the beta coefficient be assumed to be known or should model variance be included?
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
#'  \donttest{
#'  ems2 <- emulator_from_data(GillespieSIR, output_names = out_vars,
#'   ranges = ranges, c_lengths = c(0.55, 0.6, 0.59),
#'   deltas = c(0.1, 0.2, 0.2), quadratic = TRUE)
#'  ems2 # Broadly the same, but with the correlation structure modified.
#'
#'  ems2_beta <- emulator_from_data(GillespieSIR, output_names = out_vars,
#'   ranges = ranges, c_lengths = c(0.55, 0.6, 0.59),
#'   deltas = c(0.1, 0.2, 0.2), quadratic = TRUE, beta.var = TRUE)
#'  }
emulator_from_data <- function(input_data, output_names, ranges,
                               input_names = names(ranges), beta, u,
                               c_lengths, funcs, bucov, deltas, ev,
                               quadratic = TRUE, beta.var = FALSE) {
  if (missing(ranges)) {
    warning("No ranges provided; inputs assumed to be in the range [-1,1].")
    ranges <- purrr::map(input_names, ~c(-1,1))
  }
  data <- cbind(eval_funcs(scale_input, input_data[,names(ranges)], ranges), input_data[,output_names])
  if (missing(beta) || missing(funcs)) {
    if (quadratic) {
      does_add <- (choose(length(input_names)+2, length(input_names)) > nrow(input_data))
      models <- purrr::map(output_names, ~get_quadratic_model(input_data, ranges, ., add = does_add))
    }
    else {
      does_add <- (length(input_names)+1 > nrow(input_data))
      models <- purrr::map(output_names, ~get_linear_model(input_data, ranges, ., add = does_add))
    }
    all_funcs <- c(function(x) 1, purrr::map(seq_along(input_names), ~function(x) x[[.]]))
    all_coeffs <- c("(Intercept)", input_names)
    if (quadratic) {
      all_funcs <- c(all_funcs, apply(expand.grid(1:length(input_names), 1:length(input_names)), 1, function(y) function(x) x[[y[[1]]]]*x[[y[[2]]]]))
      all_coeffs <- c(all_coeffs, apply(expand.grid(input_names, input_names), 1, paste, collapse = ":"))
      all_coeffs <- sub("(.*):(\\1)", "I(\\1^2)", all_coeffs)
    }
    model_beta_mus <- purrr::map(models, ~c(.$coefficients, use.names = FALSE))
    model_basis_funcs <- purrr::map(models, ~all_funcs[all_coeffs %in% variable.names(.)])
    if (beta.var) model_beta_sigmas <- purrr::map(models, ~vcov(.))
    else model_beta_sigmas <- purrr::map(model_beta_mus, ~diag(0, nrow = length(.)))
  }
  else {
    if (missing(beta) || is.null(beta[[1]]$mu) || is.null(beta[[1]]$sigma)) stop("Basis functions provided but no regression coefficients.")
    if (missing(funcs)) stop("Regression coefficients provided but no basis functions.")
    model_beta_mus <- purrr::map(beta, ~.$mu)
    model_beta_sigmas <- purrr::map(beta, ~.$sigma)
    model_basis_funcs <- purrr::map(beta, ~funcs)
  }
  if (missing(deltas)) model_deltas <- rep(0.1, length(output_names))
  else model_deltas <- deltas
  if (missing(u) || is.null(u[[1]]$sigma) || is.null(u[[1]]$mu) || is.null(u[[1]]$corr)) {
    residuals <- purrr::map(seq_along(output_names), ~data[,output_names[[.]]] - apply(sweep(t(eval_funcs(model_basis_funcs[[.]], data[,names(ranges)])), 2, model_beta_mus[[.]], "*"), 1, sum))
    ifelse(quadratic, theta_range <- c(0.2, 1), theta_range <- c(0.2, 2))
    specs <- purrr::map(seq_along(residuals), ~get_likelihood(data[,input_names], residuals[[.]], model_basis_funcs[[.]], theta_range, delta = model_deltas[[.]]))
    model_u_sigmas <- map(specs, ~as.numeric(.$sigma))
    model_u_mus <- purrr::map(output_names, ~function(x) 0)
    if (missing(c_lengths)) c_lengths <- purrr::map(specs, ~as.numeric(.$theta))
    model_u_corr_funcs <- purrr::map(seq_along(output_names), ~function(x, xp) exp_sq(x, xp, c_lengths[[.]]))
  }
  else {
    model_u_sigmas <- purrr::map(u, ~.$sigma)
    model_u_mus <- purrr::map(u, ~.$mu)
    model_u_corr_funcs <- purrr::map(u, ~function(x,xp) .$corr(x,xp,.x$theta))
  }
  model_us <- purrr::map(seq_along(model_u_sigmas), ~list(mu = model_u_mus[[.]], sigma = model_u_sigmas[[.]], corr = model_u_corr_funcs[[.]]))
  model_betas <- purrr::map(seq_along(model_beta_mus), ~list(mu = model_beta_mus[[.]], sigma = model_beta_sigmas[[.]]))
  if (missing(bucov)) bucov <- NULL
  out_ems <- purrr::map(seq_along(model_betas), ~Emulator$new(basis_f = model_basis_funcs[[.]], beta = model_betas[[.]], u = model_us[[.]], ranges = ranges, delta = model_deltas[[.]], model = models[[.]]))
  if (missing(deltas) && !missing(ev)) {
    new_deltas <- ev/purrr::map_dbl(out_ems, ~.$u_sigma)
    new_deltas <- purrr::map_dbl(new_deltas, ~min(1/3, .))
    return(emulator_from_data(input_data, output_names, ranges, input_names, beta, u, c_lengths, funcs, bucov, new_deltas, quadratic = quadratic, beta.var = beta.var))
  }
  return(setNames(out_ems, output_names))
}

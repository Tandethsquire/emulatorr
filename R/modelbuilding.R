#' Model Generation
#'
#' Creates a best fit of coefficients for a given data set.
#'
#' There are two ways to generate the model; either start with all possible terms
#' (including cross-terms) up to order \code{n}, and then stepwise remove them;
#' or start with an intercept and stepwise add terms up to order \code{n}, only
#' retaining a term if the Information Criterion is improved. Which method is
#' chosen is dependent on the value of \code{add}; in the event where
#' \code{add = FALSE} and there are not enough degrees of freedom to start with
#' all possible terms, a warning will be given.
#'
#' @importFrom stats lm step setNames as.formula
#'
#' @param data A \code{data.frame} containing the input and output values
#' @param ranges A named \code{list} consisting of the ranges of the input parameters
#' @param output_name A string corresponding to the output to be modelled
#' @param add Should we perform stepwise add or stepwise delete? Default: \code{FALSE}
#' @param order To what order terms should the model be fitted? Default: 2 (quadratic)
#' @param u_form An upper form for the model fit. Default \code{NULL}; used internally.
#'
#' @return The fitted model
#' @export
get_coefficient_model <- function(data, ranges, output_name, add = FALSE, order = 2, u_form = NULL) {
  lower_form <- as.formula(paste(output_name, "1", sep = " ~ "))
  if (is.null(u_form)) {
    if (order == 1) {
      upper_form <- as.formula(
        paste(
          output_name, " ~ ",
          paste0(names(ranges), collapse = "+"),
          sep = ""
        )
      )
    }
    else {
      upper_form <- as.formula(
        paste(
          output_name, "~",
          paste(
            paste0(names(ranges), collapse = "+"),
            paste0("I(", names(ranges), "^2)", collapse = "+"),
            sep = "+"),
          sep = ""
        )
      )
      start_model <- get_coefficient_model(data = data, ranges = ranges, output_name = output_name, add = add, order = order, u_form = upper_form)
      a_vars <- names(start_model$coefficients)[-1]
      in_model <- purrr::map_lgl(names(ranges), ~any(grepl(., a_vars)))
      a_vars <- names(ranges)[in_model]
      upper_form <- as.formula(
        paste(
          output_name, " ~ ",
          paste(
            paste("(", paste(a_vars, collapse = "+"), ")^", order, sep = ""),
            paste0("I(", a_vars, paste(")^", order, sep = ""), collapse = "+"),
            sep = "+"
          ),
          sep = ""
        )
      )
    }
  }
  else {
    upper_form <- u_form
  }
  if (!add & (choose(length(ranges) + order, length(ranges)) > nrow(data))) {
    warning("Maximum number of regression terms is greater than the degrees of freedom in the data. Changing to add = TRUE.")
    add <- TRUE
  }
  scaled_input_data <- scale_input(data[, names(ranges)], ranges)
  full_scaled_data <- setNames(cbind(scaled_input_data, data[,output_name]), c(names(ranges), output_name))
  if (add) {
    model <- step(lm(formula = lower_form, data = full_scaled_data),
                  scope = list(lower = lower_form, upper = upper_form),
                  direction = "forward", trace = 0, k = log(nrow(data)))
  }
  else {
    model <- step(lm(formula = upper_form, data = full_scaled_data),
                  scope = list(lower = lower_form, upper = upper_form),
                  direction = "backward", trace = 0, k = log(nrow(data)))
  }
  return(model)
}

#' Hyperparameter Estimation
#'
#' Does hyperparameter fitting using the nlme package
#'
#' @importFrom nlme gls corGaus
#' @importFrom stats coef formula
#'
#' @param inputs The data to fit to
#' @param model The pre-fitted model (determined using \code{step})
#'
#' @keywords internal
#' @noRd
#'
#' @return A list of hyperparameter values.
hyperparam_fit <- function(inputs, model) {
  model_form <- formula(model$terms)
  gls_model <- gls(data = inputs, model = model_form, correlation = corGaus(nugget = TRUE))
  beta_coeffs <- c(coef(gls_model), use.names = FALSE)
  beta_var <- gls_model$varBeta
  u_sigma <- gls_model$sigma
  theta <- coef(gls_model$modelStruct)[[1]]
  delta <- coef(gls_model$modelStruct)[[2]]
  return(list(beta = list(mu = beta_coeffs, var = beta_var), sigma = u_sigma, theta = theta, delta = delta))
}

# Performs maximum likelihood estimation of the hyperparameters sigma and theta.
# At present, only does single theta estimates (i.e. the same correlation length
# is applied in all parameter directions)
get_likelihood <- function(inputs, outputs, h, theta_range, n_ints = 100, delta = 0.1) {
  best_L <- -Inf
  best_sigma <- best_theta <- NULL
  hs <- apply(inputs, 1, function(y) purrr::map_dbl(h, purrr::exec, y))
  corr_func <- function(x, xp, theta, delta) {
    (1-delta) * exp(-sum((x-xp)^2/theta^2)) + ifelse(all(x == xp), delta, 0)
  }
  corr_matrix <- function(points, theta, delta) {
    apply(points, 1, function(x) apply(points, 1, function(y) corr_func(x,y,theta,delta)))
  }
  t_list <- seq(theta_range[1], theta_range[2], length.out = n_ints)
  likelihood_list <- lapply(t_list, function(i) {
    cmat <- corr_matrix(inputs, i, delta)
    cm <- minv(cmat)
    #beta_est <- chol2inv(chol(hs %*% cm %*% t(hs))) %*% hs %*% cm %*% outputs
    beta_est <- rep(0, nrow(hs))
    var_est <- t(outputs - t(hs) %*% beta_est) %*% cm %*% (outputs - t(hs) %*% beta_est)
    lik <- det(cmat)^(-1/2) * det(mmult(hs, mmult(cm, t(hs))))^(-1/2) * var_est^(-(nrow(inputs)-nrow(hs))/2)
    return(list(L = lik, S = sqrt(var_est), Th = i))
  })
  liks <- purrr::map_dbl(likelihood_list, ~.$L)
  l_ind <- which.max(liks)
  return(list(sigma = likelihood_list[[l_ind]]$S/sqrt(length(outputs)), theta = likelihood_list[[l_ind]]$Th))
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
#'   ranges = ranges, quadratic = FALSE)
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
      models <- purrr::map(output_names, ~get_coefficient_model(input_data, ranges, ., add = does_add))
    }
    else {
      does_add <- (length(input_names)+1 > nrow(input_data))
      models <- purrr::map(output_names, ~get_coefficient_model(input_data, ranges, ., add = does_add, order = 1))
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
    estimate_method <- 'nl'
    if (estimate_method == "my") {
      residuals <- purrr::map(models, ~.$residuals)
      ifelse(quadratic, theta_range <- c(0.2, 1), theta_range <- c(0.2, 2))
      specs <- purrr::map(seq_along(residuals), ~get_likelihood(data[,input_names], residuals[[.]], model_basis_funcs[[.]], theta_range, delta = model_deltas[[.]]))
    }
    else {
      specs <- purrr::map(seq_along(output_names), ~hyperparam_fit(data[,c(input_names, output_names[[.]])], models[[.]]))
      model_u_sigmas <- purrr::map(specs, ~as.numeric(.$sigma))
      model_u_mus <- purrr::map(output_names, ~function(x) 0)
      deltas <- purrr::map_dbl(specs, ~.$delta)
    }
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
  for (i in 1:length(out_ems)) out_ems[[i]]$output_name <- output_names[[i]]
  names(out_ems) <- output_names
  return(out_ems)
}

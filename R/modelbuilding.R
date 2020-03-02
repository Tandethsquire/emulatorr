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
  scaled_input_data <- t(apply(data[,names(ranges)], 1, scale_input, ranges))
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
  scaled_input_data <- t(apply(data[,names(ranges)], 1, scale_input, ranges))
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

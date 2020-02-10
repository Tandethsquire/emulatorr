Correlator <- R6::R6Class(
  "Correlator",
  public = list(
    cov_func = NULL,
    exp_func = NULL,
    initialize = function(cov, mu) {
      if (length(formals(cov)) != 2) stop("Covariance function must take exactly two arguments.")
      if (length(formals(mu)) != 1) stop("Expectation function must take exactly one argument.")
      if (length(cov(1, 1)) != 1 | length(mu(1)) != 1) stop("The covariance and expectation functions must return a single value.")
      self$cov_func = cov
      self$exp_func = mu
    },
    get_exp =  function(x) {
      return(self$exp_func(x))
    },
    get_cov = function(x, xp = NULL) {
      return(ifelse(is.null(xp), self$cov_func(x, x), self$cov_func(x, xp)))
    }
  )
)

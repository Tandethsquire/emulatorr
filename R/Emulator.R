Emulator <- R6::R6Class(
  "Emulator",
  public = list(
    beta = NULL,
    u = NULL,
    beta_u_cov = NULL,
    basis_f = NULL,
    param_ranges = NULL,
    initialize = function(funcs, beta, u, bucov=NULL, n_inputs=NULL, ranges=purrr::map(seq_along(1:n_inputs), ~c(-1,1))) {
      if (class(u)[1] != "Correlator") stop("u must be a correlator object.")
      self$basis_f = funcs
      self$beta = beta
      self$u = u
      self$beta_u_cov = ifelse(is.null(bucov), function(x) rep(0, length(beta$mu)), bucov)
      self$param_ranges = ranges
    },
    get_exp = function(x) {
      x <- private$scale_inputs(x)
      purrr::map_dbl(self$basis_f, purrr::exec, x) %*% self$beta$mu  +self$u$get_exp(x)
    },
    get_var = function(x) {
      x <- private$scale_inputs(x)
      f_map <- purrr::map_dbl(self$basis_f, purrr::exec, x)
      return(f_map %*% self$beta$sigma %*% f_map
        + 2 * f_map %*% self$beta_u_cov(x)
        + self$u$get_cov(x))
    },
    bayes_adjust = function(inputs, outputs) {
      inputs <- t(apply(inputs, 1, private$scale_inputs))
      omega <- apply(inputs, 1, function(x) apply(inputs, 1, function(y) self$u$get_cov(x,y)))
      O <- solve(omega)
      G <- t(apply(inputs, 1, function(x) purrr::map_dbl(self$basis_f, ~purrr::exec(.x, c(x)))))
      GOG <- t(G) %*% O %*% G
      gls <- solve(GOG) %*% t(G) %*% O %*% outputs
      ifelse(all(self$beta$sigma == 0), siginv <- self$beta$sigma, siginv <- solve(self$beta$sigma))
      new_beta_var <- solve(GOG + siginv)
      new_beta_exp <- c(new_beta_var %*% (GOG %*% gls + siginv %*% self$beta$mu))
      point_cov <- function(x) apply(inputs, 1, function(y) self$u$get_cov(x, y)/sqrt(self$u$get_cov(x) * self$u$get_cov(y)))
      new_u_exp <- function(x) {
        return(self$u$get_exp(x)
               + (self$beta_u_cov(x) %*% t(G) + self$u$get_cov(x) * point_cov(x)) %*% O %*% (outputs - G %*% new_beta_exp))
      }
      new_u_var <- function(x, y) {
        correlations <- function(x) G %*% self$beta_u_cov(x) + self$u$get_cov(x) * point_cov(x)
        return(self$u$get_cov(x)
               -t(correlations(x)) %*% (O - O %*% G %*% new_beta_var %*% t(G) %*% O) %*% correlations(y))
      }
      new_beta_u_cov <- function(x) {
        siginv %*% new_beta_var %*% self$beta_u_cov(x) - new_beta_var %*% t(G) %*% O %*% point_cov(x) * self$u$get_cov(x)
      }
      new_u <- Correlator$new(new_u_var, new_u_exp)
      return(Emulator$new(self$basis_f, list(mu = new_beta_exp, sigma = new_beta_var), new_u, new_beta_u_cov, n_inputs = self$n_inputs, ranges = self$param_ranges))
    },
    implausibility = function(x, z) {
      x <- private$scale_inputs(x)
      if (is.numeric(z))
        output <- list(val = z, sigma = 0)
      else output <- z
      imp_var <- self$get_var(x) + output$sigma^2
      return(sqrt((output$val-self$get_exp(x))^2/imp_var))
    },
    print = function(...) {
      cat("Emulator: \n")
      cat("Beta: \n")
      cat("\t Mu: ", paste(self$beta$mu, collapse="; "), "\n")
      cat("\t Sigma:", "matrix(", paste(self$beta$sigma, collapse=", "), ")\n")
    }
  ),
  private = list(
    scale_inputs = function(x, scale = TRUE) {
      centers <- purrr::map_dbl(self$param_ranges, ~(.x[2]+.x[1])/2)
      scales <- purrr::map_dbl(self$param_ranges, ~(.x[2]-.x[1])/2)
      if (scale) return((x - centers)/scales)
      return(x * scales + centers)
    }
  )
)

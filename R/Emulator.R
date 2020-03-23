Emulator <- R6::R6Class(
  "Emulator",
  private = list(
    data_corrs = NULL,
    design_matrix = NULL,
    u_exp_modifier = NULL,
    u_var_modifier = NULL,
    beta_u_cov_modifier = NULL
  ),
  public = list(
    basis_f = NULL,
    beta_mu = NULL,
    beta_sigma = NULL,
    u_mu = NULL,
    u_sigma = NULL,
    corr_func = NULL,
    beta_u_cov = NULL,
    in_data = NULL,
    out_data = NULL,
    ranges = NULL,
    delta = NULL,
    initialize = function(basis_f, beta, u, ranges, bucov = NULL, data = NULL, delta = 0) {
      self$basis_f <- basis_f
      self$beta_mu <- beta$mu
      self$beta_sigma <- beta$sigma
      self$u_mu <- u$mu
      self$u_sigma <- u$sigma
      self$delta <- delta
      self$corr_func <- function(x, xp) {
        if (sum((x-xp)^2) < 0.0001) extra <- 1
        else extra <- 0
        (1-delta) * u$corr(x, xp) + delta * extra
      }
      if (is.null(bucov))
        self$beta_u_cov <- function(x) rep(0, length(self$beta_mu))
      else
        self$beta_u_cov <- bucov
      self$ranges <- ranges
      if (is.null(ranges)) stop("Ranges for the parameters must be specified.")
      if (!is.null(data)) {
        self$in_data <- eval_funcs(scale_input, data[,names(self$ranges)], self$ranges)
        self$out_data <- data[, !(names(data) %in% names(self$ranges))]
      }
      if (!is.null(self$in_data))
      {
        private$data_corrs <- chol2inv(chol(self$u_sigma^2 * apply(self$in_data, 1, function(y) eval_funcs(self$corr_func, self$in_data, y))))
        private$design_matrix <- t(eval_funcs(self$basis_f, self$in_data))
        private$u_var_modifier <- private$data_corrs %*% private$design_matrix %*% self$beta_sigma %*% t(private$design_matrix) %*% private$data_corrs
        private$u_exp_modifier <- private$data_corrs %*% (self$out_data - private$design_matrix %*% self$beta_mu)
        private$beta_u_cov_modifier <- self$beta_sigma %*% t(private$design_matrix) %*% private$data_corrs
      }
    },
    get_exp = function(x) {
      x <- eval_funcs(scale_input, x, self$ranges)
      g <- t(eval_funcs(self$basis_f, x))
      bu <- eval_funcs(self$beta_u_cov, x)
      beta_part <- g %*% self$beta_mu
      u_part <- eval_funcs(self$u_mu, x)
      if (!is.null(self$in_data)) {
        c_data <- apply(self$in_data, 1, function(y) eval_funcs(self$corr_func, x, y))
        u_part <- u_part + (t(bu) %*% t(private$design_matrix) + self$u_sigma^2 * c_data) %*% private$u_exp_modifier
      }
      return(beta_part + u_part)
    },
    get_cov = function(x, xp = NULL, full = FALSE) {
      x <- eval_funcs(scale_input, x, self$ranges)
      if (is.null(xp)) xp <- x
      else xp <- eval_funcs(scale_input, xp, ranges)
      if (full) {
        if (is.null(names(xp)) || length(xp[,names(xp)[1]]) == 1) x_xp_c <- eval_funcs(self$corr_func, x, xp)
        else x_xp_c <- apply(xp, 1, function(y) eval_funcs(self$corr_func, x, y))
        g_x <- t(eval_funcs(self$basis_f, x))
        g_xp <- t(eval_funcs(self$basis_f, xp))
        beta_part <- g_x %*% self$beta_sigma %*% t(g_xp)
        u_part <- self$u_sigma^2 * x_xp_c
        bupart_x <- eval_funcs(self$beta_u_cov, x)
        bupart_xp <- eval_funcs(self$beta_u_cov, xp)
        if (!is.null(self$in_data)) {
          if (is.null(names(x)) || length(x[,names(x)[1]]) == 1) c_x <- t(eval_funcs(self$corr_func, self$in_data, x))
          else c_x <- apply(self$in_data, 1, function(y) eval_funcs(self$corr_func, x, y))
          if (is.null(names(xp)) || length(xp[,names(xp)[1]]) == 1) c_xp <- t(eval_funcs(self$corr_func, self$in_data, xp))
          else c_xp <- apply(self$in_data, 1, function(y) eval_funcs(self$corr_func, xp, y))
          u_part <- u_part - self$u_sigma^2 * c_x %*% (private$data_corrs - private$u_var_modifier) %*% t(c_xp) * self$u_sigma^2
          bupart_x <- bupart_x - private$beta_u_cov_modifier %*% (private$design_matrix %*% bupart_x + self$u_sigma^2 * t(c_x))
          bupart_xp <- bupart_xp - private$beta_u_cov_modifier %*% (private$design_matrix %*% bupart_xp + self$u_sigma * t(c_xp))
        }
        bupart <- g_x %*% bupart_xp + t(bupart_x) %*% t(g_xp)
      }
      else {
        if (length(x) != length(xp)) stop("Can't compute the diagonal elements of a non-square covariance matrix.")
        if (class(x) == 'numeric' || length(x[,names(x)[1]]) == 1) {
          beta_part <- t(eval_funcs(self$basis_f, x)) %*% self$beta_sigma %*% eval_funcs(self$basis_f, xp)
          u_part <- self$u_sigma^2 * self$corr_func(x, xp)
          bupart_x <- eval_funcs(self$beta_u_cov, x)
          bupart_xp <- eval_funcs(self$beta_u_cov, xp)
          if (!is.null(self$in_data)) {
            c_x <- t(eval_funcs(self$corr_func, self$in_data, x))
            c_xp <- t(eval_funcs(self$corr_func, self$in_data, xp))
            u_part <- u_part - self$u_sigma^2 * c_x %*% (private$data_corrs - private$u_var_modifier) %*% t(c_xp) * self$u_sigma^2
            bupart_x <- bupart_x - private$beta_u_cov_modifier %*% (private$design_matrix %*% bupart_x + self$u_sigma^2 *t(c_x))
            bupart_xp <- bupart_xp - private$beta_u_cov_modifier %*% (private$design_matrix %*% bupart_xp + self$u_sigma^2 *t(c_xp))
          }
          bupart <- t(eval_funcs(self$basis_f, x)) %*% bupart_xp + t(bupart_x) %*% eval_funcs(self$basis_f, xp)
        }
        else {
          xseq <- seq_along(x[,1])
          beta_part <- purrr::map_dbl(xseq, ~t(eval_funcs(self$basis_f, x[.x,])) %*% self$beta_sigma %*% eval_funcs(self$basis_f, xp[.x,]))
          u_part <- purrr::map_dbl(xseq, ~self$corr_func(x[.x,], xp[.x,])) * self$u_sigma^2
          bupart_x <- eval_funcs(self$beta_u_cov, x)
          bupart_xp <- eval_funcs(self$beta_u_cov, xp)
          if (!is.null(self$in_data)) {
            c_x <- t(apply(x, 1, function(y) eval_funcs(self$corr_func, self$in_data, y)))
            c_xp <- t(apply(xp, 1, function(y) eval_funcs(self$corr_func, self$in_data, y)))
            u_part <- u_part - self$u_sigma^2 * rowSums((c_x %*% (private$data_corrs - private$u_var_modifier)) * c_xp) * self$u_sigma^2
            bupart_x <- bupart_x - private$beta_u_cov_modifier %*% (private$design_matrix %*% bupart_x + self$u_sigma^2 * t(c_x))
            bupart_xp <- bupart_xp - private$beta_u_cov_modifier %*% (private$design_matrix %*% bupart_xp + self$u_sigma^2 * t(c_xp))
          }
          bupart <- purrr::map_dbl(xseq, ~t(eval_funcs(self$basis_f, x[.x,])) %*% bupart_xp[,.x] + t(bupart_x[,.x]) %*% eval_funcs(self$basis_f, xp[.x,]))
        }
      }
      return(beta_part + u_part + bupart)
    },
    implausibility = function(x, z) {
      if (is.numeric(z)) output <- list(val = z, sigma = 0)
      else output <- z
      imp_var <- self$get_cov(x) + output$sigma^2
      return(sqrt((output$val - self$get_exp(x))^2/imp_var))
    },
    adjust = function(data, out_name) {
      this_data_in <- eval_funcs(scale_input, data[,names(self$ranges)], self$ranges)
      this_data_out <- data[,out_name]
      G <- t(eval_funcs(self$basis_f, this_data_in))
      O <- chol2inv(chol(self$u_sigma^2 * apply(this_data_in, 1, function(y) eval_funcs(self$corr_func, this_data_in, y))))
      if (all(eigen(self$beta_sigma)$values == 0))
      {
        new_beta_var = self$beta_sigma
        new_beta_exp = self$beta_mu
      }
      else
      {
        siginv <- chol2inv(chol(self$beta_sigma))
        new_beta_var <- chol2inv(chol(t(G) %*% O %*% G + siginv))
        new_beta_exp <- new_beta_var %*% (siginv %*% self$beta_mu + t(G) %*% O %*% this_data_out)
      }
      new_em <- Emulator$new(self$basis_f, list(mu = new_beta_exp, sigma = new_beta_var),
                             u = list(mu = self$u_mu, sigma = self$u_sigma, corr = self$corr_func),
                             bucov = self$beta_u_cov, ranges = self$ranges, data = data[,c(names(self$ranges), out_name)], delta = self$delta)
      return(new_em)
    },
    print = function(...) {
      cat("Parameters and ranges: ", paste(names(self$ranges), paste0(self$ranges), sep = ": ", collapse= "; "), "\n")
      cat("Specifications: \n")
      cat("\t Basis functions: ", paste(purrr::map(self$basis_f, ~function_to_names(.x, names(self$ranges))), collapse = "; "), "\n")
      cat("\t Beta Expectation: ", paste(round(self$beta_mu, 4), collapse = "; "), "\n")
      cat("\t Beta Variance (eigenvalues): ", paste(round(eigen(self$beta_sigma)$values, 4), collapse = "; "), "\n")
      cat("Correlation Structure: \n")
      if (!is.null(private$data_corrs)) cat("Non-stationary covariance - prior specifications below \n")
      cat("\t Variance: ", self$u_sigma^2, "\n")
      cat("\t Expectation: ", self$u_mu(rep(0, length(ranges))), "\n")
      cat("\t Correlation length: ", sqrt(-0.25/(log(self$corr_func(0,0.5))-log(1-self$delta))), "\n")
      cat("\t Nugget term: ", self$delta, "\n")
      cat("Mixed covariance: ", self$beta_u_cov(rep(0, length(ranges))), "\n")
    }
  )
)

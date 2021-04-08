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
    theta = NULL,
    ranges = NULL,
    delta = NULL,
    model = NULL,
    model_terms = NULL,
    active_vars = NULL,
    corr = NULL,
    o_em = NULL,
    output_name = NULL,
    data_diag = NULL,
    initialize = function(basis_f, beta, u, ranges, bucov = NULL, data = NULL, delta = 0, model = NULL, original_em = NULL, out_name = NULL, a_vars = NULL, data_diag = 0) {
      self$model <- model
      self$model_terms <- tryCatch(
        c("1", labels(terms(self$model))),
        error = function(e) return(NULL)
      )
      self$o_em <- original_em
      self$basis_f <- basis_f
      self$beta_mu <- beta$mu
      self$beta_sigma <- beta$sigma
      self$u_mu <- u$mu
      self$u_sigma <- u$sigma
      self$delta <- delta
      self$corr <- u$corr
      self$data_diag <- data_diag
      if (!is.null(out_name)) self$output_name <- out_name
      if (is.null(a_vars)) {
        self$active_vars <- purrr::map_lgl(seq_along(ranges), function(x) {
          point_vec <- c(rep(0, x-1), 1, rep(0, length(ranges)-x))
          func_vals <- purrr::map_dbl(self$basis_f, purrr::exec, point_vec)
          sum(func_vals) > 1
        })
      }
      else self$active_vars <- a_vars
      self$corr_func <- function(x, xp) {
        extra <- if(sum((x-xp)^2) > 1e-6) 0 else 1
        (1-delta) * self$corr(x[self$active_vars], xp[self$active_vars]) + delta * extra
      }
      self$beta_u_cov <- if (is.null(bucov)) function(x) rep(0, length(self$beta_mu)) else bucov
      self$theta <- sqrt(-1/log(self$corr(0,0.5)))/2
      if (is.null(ranges)) stop("Ranges for the parameters must be specified.")
      self$ranges <- ranges
      if (!is.null(data)) {
        self$in_data <- data.matrix(eval_funcs(scale_input, data[,names(self$ranges)], self$ranges))
        self$out_data <- data[, !(names(data) %in% names(self$ranges))]
      }
      if (!is.null(self$in_data)) {
        private$data_corrs <- chol2inv(chol(self$u_sigma^2 * apply(self$in_data, 1, function(y) apply(self$in_data, 1, self$corr_func, y)) + diag(self$data_diag, nrow = nrow(self$in_data))))
        private$design_matrix <- t(apply(self$in_data, 1, function(x) purrr::map_dbl(self$basis_f, purrr::exec, x)))
        if (nrow(private$design_matrix) == 1) private$design_matrix <- t(private$design_matrix)
        private$u_var_modifier <- private$data_corrs %*% private$design_matrix %*% self$beta_sigma %*% t(private$design_matrix) %*% private$data_corrs
        private$u_exp_modifier <- private$data_corrs %*% (self$out_data - private$design_matrix %*% self$beta_mu)
        private$beta_u_cov_modifier <- self$beta_sigma %*% t(private$design_matrix) %*% private$data_corrs
      }
    },
    get_exp = function(x, p = NULL, local.var = TRUE) {
      if (is.null(self$model_terms) && !is.null(p)) {
        warning("Can't extract functional form of basis functions, so cannot differentiate them. Setting p = NULL.")
        p = NULL
      }
      x <- x[, names(x) %in% names(self$ranges)]
      x <- eval_funcs(scale_input, x, self$ranges)
      if (is.null(p)) {
        g <- t(apply(x, 1, function(y) purrr::map_dbl(self$basis_f, purrr::exec, y)))
      }
      else {
        if (!p %in% names(self$ranges)) return(0)
        g_d <- purrr::map(self$model_terms, ~D(parse(text = sub("I\\((\\w*\\^\\d*)\\)", "\\1", gsub(":", "*", .))), p))
        g <- matrix(unlist(do.call('rbind', purrr::map(seq_along(x[,1]), function(y) purrr::map(g_d, ~eval(., envir = x[y,]))))), nrow = nrow(x))
      }
      x <- data.matrix(x)
      ## Need to think about how Cov[beta, u] will change if we use derivatives.
      bu <- t(apply(x, 1, self$beta_u_cov))
      if (length(self$beta_mu) == 1) beta_part <- g * self$beta_mu
      else beta_part <- g %*% self$beta_mu
      u_part <- apply(x, 1, self$u_mu)
      if (!is.null(self$in_data)) {
        if (!is.null(p) && local.var) {
          c_data <- -2*apply(self$in_data, 1, function(y) x[,p] - y[p]) * apply(self$in_data, 1, function(y) apply(x, 1, self$corr_func, y))/self$theta^2
        }
        else {
          c_data <- apply(self$in_data, 1, function(y) apply(x, 1, self$corr_func, y))
        }
        if (length(self$beta_mu) == 1) {
          u_part <- t(u_part + (t(bu) %*% t(private$design_matrix) + self$u_sigma^2 * c_data) %*% private$u_exp_modifier)
        }
        else
          u_part <- u_part + (bu %*% t(private$design_matrix) + self$u_sigma^2 * c_data) %*% private$u_exp_modifier
      }
      if (length(self$beta_mu) == 1) return(c(beta_part + u_part))
      return(beta_part + u_part)
    },
    get_cov = function(x, p = NULL, xp = NULL, full = FALSE, pp = NULL, local.var = TRUE) {
      if (is.null(self$model_terms) && !is.null(p)) {
        warning("Can't extract functional form of basis functions, so cannot differentiate them. Setting p = NULL.")
        p = NULL
      }
      x <- x[, names(x) %in% names(self$ranges)]
      x <- eval_funcs(scale_input, x, self$ranges)
      if (!is.null(p)) {
        if (is.null(pp)) pp <- p
        if (!p %in% names(self$ranges) || !pp %in% names(self$ranges)) {
          if (!full) return(0)
          else return(matrix(0, nrow = nrow(x), ncol = nrow(xp)))
        }
        g_d <- purrr::map(self$model_terms, ~D(parse(text = sub("I\\((\\w*\\^\\d*)\\)", "\\1", gsub(":", "*", .))), p))
        gp_d <- purrr::map(self$model_terms, ~D(parse(text = sub("I\\((\\w*\\^\\d*)\\)", "\\1", gsub(":", "*", .))), pp))
        g_x <- t(matrix(unlist(do.call('rbind', purrr::map(seq_along(x[,1]), function(y) purrr::map(g_d, ~eval(., envir = x[y,]))))), nrow = nrow(x)))
        gp_x <- t(matrix(unlist(do.call('rbind', purrr::map(seq_along(x[,1]), function(y) purrr::map(gp_d, ~eval(., envir = x[y,]))))), nrow = nrow(x)))
      }
      else {
        g_x <- apply(x, 1, function(y) purrr::map_dbl(self$basis_f, purrr::exec, y))
      }
      x <- data.matrix(x)
      bupart_x <- apply(x, 1, self$beta_u_cov)
      if (is.null(xp)) {
        xp <- x
        if (!is.null(pp)) g_xp <- gp_x
        else g_xp <- g_x
        bupart_xp <- bupart_x
      }
      else {
        if (!is.null(p)) {
          xp <- eval_funcs(scale_input, xp, self$ranges)
          g_xp <- t(matrix(unlist(do.call('rbind', purrr::map(seq_along(xp[,1]), function(y) purrr::map(gp_d, ~eval(., envir = xp[y,]))))), nrow = nrow(xp)))
        }
        else {
          g_xp <- apply(xp, 1, function(y) purrr::map_dbl(self$basis_f, purrr::exec, y))
        }
        xp <- data.matrix(xp)
        bupart_xp <- apply(xp, 1, self$beta_u_cov)
      }
      if (full || nrow(x) != nrow(xp)) {
        if (!is.null(p) && local.var) {
          if (p == pp)
            x_xp_c <- (2/self$theta^2 - 4/self$theta^4 * apply(xp, 1, function(y) apply(x, 1, function(x) (x[p] - y[p])^2))) * apply(xp, 1, function(y) apply(x, 1, self$corr_func, y))
          else
            x_xp_c <- -4/self$theta^4 * apply(xp, 1, function(y) apply(x, 1, function(x) (x[p]-y[p])*(x[pp]-y[pp]))) * apply(xp, 1, function(y) apply(x, 1, self$corr_func, y))
        }
        else {
          x_xp_c <- apply(xp, 1, function(y) apply(x, 1, self$corr_func, y))
        }
        if (is.null(nrow(g_x))) beta_part <- g_x %*% self$beta_sigma %*% g_xp
        else beta_part <- t(g_x) %*% self$beta_sigma %*% g_xp
        u_part <- self$u_sigma^2 * x_xp_c
        if (!is.null(self$in_data)) {
          if (!is.null(p) && local.var) {
            c_x <- -2*apply(self$in_data, 1, function(y) x[,p] - y[p]) * apply(self$in_data, 1, function(y) apply(x, 1, self$corr_func, y))/self$theta^2
            c_xp <- -2*apply(self$in_data, 1, function(y) xp[,pp] - y[pp]) * apply(self$in_data, 1, function(y) apply(xp, 1, self$corr_func, y))/self$theta^2
          }
          else {
            c_x <- apply(self$in_data, 1, function(y) apply(x, 1, self$corr_func, y))
            c_xp <- apply(self$in_data, 1, function(y) apply(xp, 1, self$corr_func, y))
          }
          if (nrow(x) == 1) {
            c_x <- t(c_x)
            c_xp <- t(c_xp)
          }
          u_part <- u_part - self$u_sigma^4 * c_x %*% (private$data_corrs - private$u_var_modifier) %*% t(c_xp)
          bupart_x <- bupart_x - private$beta_u_cov_modifier %*% (private$design_matrix %*% bupart_x + self$u_sigma^2 * t(c_x))
          bupart_xp <- bupart_xp - private$beta_u_cov_modifier %*% (private$design_matrix %*% bupart_xp + self$u_sigma^2 * t(c_xp))
          if(is.null(nrow(g_x))) {
            bupart_x <- c(bupart_x)
            bupart_xp <- c(bupart_xp)
          }
        }
        if (is.null(nrow(g_x))) bupart <- outer(g_x, bupart_xp, "*") + outer(bupart_x, g_xp, "*")
        else bupart <- t(g_x) %*% bupart_xp + t(bupart_x) %*% g_xp
      }
      else {
        point_seq <- 1:nrow(x)
        if (is.null(nrow(g_x))) beta_part <- diag(diag(self$beta_sigma) * outer(g_x, g_xp))
        else beta_part <- purrr::map_dbl(point_seq, ~g_x[,.] %*% self$beta_sigma %*% g_xp[,.])
        if (!is.null(p) && local.var) {
          if (p == pp)
            u_part <- self$u_sigma^2 * purrr::map_dbl(point_seq, ~(2/self$theta^2 - 4/self$theta^4 * (x[.,p] - xp[.,p])^2) * self$corr_func(x[.,], xp[.,]))
          else
            u_part <- self$u_sigma^2 * purrr::map_dbl(point_seq, ~-4/self$theta^4 * (x[.,p]-xp[.,p]) * (x[.,pp]-xp[.,pp]) * self$corr_func(x[.,], xp[.,]))
        }
        else {
          u_part <- purrr::map_dbl(point_seq, ~self$corr_func(x[.,], xp[.,])) * self$u_sigma^2
        }
        if (!is.null(self$in_data)) {
          if (!is.null(p) && local.var) {
            c_x <- -2*apply(self$in_data, 1, function(y) x[,p] - y[p]) * apply(self$in_data, 1, function(y) apply(x, 1, self$corr_func, y))/self$theta^2
            c_xp <- -2*apply(self$in_data, 1, function(y) xp[,pp] - y[pp]) * apply(self$in_data, 1, function(y) apply(xp, 1, self$corr_func, y))/self$theta^2
          }
          else {
            c_x <- apply(self$in_data, 1, function(y) apply(x, 1, self$corr_func, y))
            c_xp <- if(all(x == xp)) c_x else apply(self$in_data, 1, function(y) apply(xp, 1, self$corr_func, y))
          }
          if (nrow(x) == 1) {
            c_x <- t(c_x)
            c_xp <- t(c_xp)
          }
          u_part <- u_part - self$u_sigma^4 * rowSums((c_x %*% (private$data_corrs - private$u_var_modifier)) * c_xp)
          bupart_x <- bupart_x - private$beta_u_cov_modifier %*% (private$design_matrix %*% bupart_x + self$u_sigma^2 * t(c_x))
          bupart_xp <- bupart_xp - private$beta_u_cov_modifier %*% (private$design_matrix %*% bupart_xp + self$u_sigma^2 * t(c_xp))
        }
        if (is.null(nrow(g_x))) bupart <- purrr::map_dbl(point_seq, ~t(purrr::map_dbl(self$basis_f, purrr::exec, x[.,])) * bupart_xp[.] + t(bupart_x[.] * purrr::map_dbl(self$basis_f, purrr::exec, xp[.,])))
        else bupart <- purrr::map_dbl(point_seq, ~t(purrr::map_dbl(self$basis_f, purrr::exec, x[.,])) %*% bupart_xp[,.x] + t(bupart_x[,.x] %*% purrr::map_dbl(self$basis_f, purrr::exec, xp[.,])))
      }
      final_out <- beta_part + u_part + bupart
      return(round(final_out,10))
    },
    implausibility = function(x, z, cutoff = NULL) {
      output <- if (is.numeric(z)) list(val = z, sigma = 0) else z
      imp_var <- self$get_cov(x) + output$sigma^2
      imp <- sqrt((output$val - self$get_exp(x))^2/imp_var)
      if (is.null(cutoff)) return(imp)
      return(imp<=cutoff)
    },
    adjust = function(data, out_name) {
      this_data_in <- data.matrix(eval_funcs(scale_input, data[,names(self$ranges)], self$ranges))
      this_data_out <- data[,out_name]
      if (all(eigen(self$beta_sigma)$values == 0)) {
        new_beta_var <- self$beta_sigma
        new_beta_exp <- self$beta_mu
      }
      else {
        G <- apply(this_data_in, 1, function(x) purrr::map_dbl(self$basis_f, purrr::exec, x))
        O <- chol2inv(chol(apply(this_data_in, 1, function(x) apply(this_data_in, 1, self$corr_func, x)) + diag(self$data_diag, nrow = nrow(this_data_in))))
        siginv <- chol2inv(chol(self$beta_sigma))
        new_beta_var <- chol2inv(chol(G %*% O %*% (if(is.null(nrow(G))) G else t(G)) + siginv))
        new_beta_exp <- new_beta_var %*% (siginv %*% self$beta_mu + G %*% O %*% this_data_out)
      }
      new_em <- Emulator$new(self$basis_f, list(mu = new_beta_exp, sigma = new_beta_var),
                             u = list(mu = self$u_mu, sigma = self$u_sigma, corr = self$corr),
                             bucov = self$beta_u_cov, ranges = self$ranges, data = data[, c(names(self$ranges),out_name)],
                             delta = self$delta, original_em = self, out_name = out_name, model = self$model,
                             a_vars = self$active_vars, data_diag = self$data_diag)
      return(new_em)
    },
    get_hessian = function(x, local.var = TRUE) {
      x <- x[, names(x) %in% names(self$ranges)]
      x <- eval_funcs(scale_input, x, self$ranges)
      result <- do.call('rbind', purrr::map(names(self$ranges), function(p1) {
        purrr::map_dbl(names(self$ranges), function(p2) {
          g_d <- purrr::map(self$model_terms, ~D(D(parse(text = sub("I\\((\\w*\\^\\d*)\\)", "\\1", gsub(":", "*", .))), p1),p2))
          g <- matrix(unlist(do.call('rbind', purrr::map(seq_along(x[,1]), function(y) purrr::map(g_d, ~eval(., envir = x[y,]))))), nrow = nrow(x))
          x <- data.matrix(x)
          ## Need to think about how Cov[beta, u] will change if we use derivatives.
          bu <- t(apply(x, 1, self$beta_u_cov))
          beta_part <- g %*% self$beta_mu
          u_part <- apply(x, 1, self$u_mu)
          if (!is.null(self$in_data)) {
            if (local.var) {
              if (p1 == p2)
                c_data <- (4/self$theta^4 * apply(self$in_data, 1, function(y) (x[,p1]-y[p1])^2) - 2/self$theta^2) * apply(self$in_data, 1, function(y) apply(x, 1, self$corr_func, y))
              else
                c_data <- 4/self$theta^4 * apply(self$in_data, 1, function(y) (x[,p1]-y[p1])*(x[,p2]-y[p2])) * apply(self$in_data, 1, function(y) apply(x, 1, self$corr_func, y))
            }
            else {
              c_data <- apply(self$in_data, 1, function(y) apply(x, 1, self$corr_func, y))
            }
            u_part <- u_part + (bu %*% t(private$design_matrix) + self$u_sigma^2 * c_data) %*% private$u_exp_modifier
          }
          return(beta_part + u_part)
        })
      }))
      return(result)
    },
    set_sigma = function(sigma) {
      if(is.null(self$o_em)) {
        new_em <- self$clone()
        new_em$u_sigma <- sigma
        return(new_em)
      }
      new_o_em <- self$o_em$clone()
      new_o_em$u_sigma <- sigma
      dat <- setNames(data.frame(cbind(eval_funcs(scale_input, data.frame(self$in_data), self$ranges, FALSE), self$out_data)), c(names(self$ranges), self$output_name))
      return(new_o_em$adjust(dat, self$output_name))
    },
    ## This will not work properly if the correlation structure isn't exp_sq
    set_theta = function(theta) {
      if(is.null(self$o_em)) {
        new_em <- self$clone()
        new_em$corr <- function(x, xp) exp_sq(x, xp, theta)
        return(new_em)
      }
      new_o_em <- self$o_em$clone()
      new_o_em$corr <- function(x, xp) exp_sq(x, xp, theta)
      dat <- setNames(data.frame(cbind(eval_funcs(scale_input, data.frame(self$in_data), self$ranges, FALSE), self$out_data)), c(names(self$ranges), self$output_name))
      return(new_o_em$adjust(dat, self$output_name))
    },
    print = function(...) {
      cat("Parameters and ranges: ", paste(names(self$ranges), paste0(purrr::map(self$ranges, round, 4)), sep = ": ", collapse= "; "), "\n")
      cat("Specifications: \n")
      if (is.null(self$model))
        cat("\t Basis functions: ", paste0(names(self$o_em$model$coefficients), collapse="; "), "\n")
      else
        cat("\t Basis Functions: ", paste0(names(self$model$coefficients), collapse="; "), "\n")
      cat("\t Active variables: ", paste0(names(self$ranges)[self$active_vars], collapse="; "), "\n")
      cat("\t Beta Expectation: ", paste(round(self$beta_mu, 4), collapse = "; "), "\n")
      cat("\t Beta Variance (eigenvalues): ", paste(round(eigen(self$beta_sigma)$values, 4), collapse = "; "), "\n")
      cat("Correlation Structure: \n")
      if (!is.null(private$data_corrs)) cat("Non-stationary covariance - prior specifications below \n")
      cat("\t Variance: ", self$u_sigma^2, "\n")
      cat("\t Expectation: ", self$u_mu(rep(0, length(ranges))), "\n")
      cat("\t Correlation length: ", self$theta, "\n")
      cat("\t Nugget term: ", self$delta, "\n")
      cat("Mixed covariance: ", self$beta_u_cov(rep(0, length(ranges))), "\n")
    }
  )
)

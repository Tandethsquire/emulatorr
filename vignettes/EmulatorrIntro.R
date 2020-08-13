## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 12,
  fig.height= 8
)
getOutputs <- function(points, times, ev = TRUE) {
  out_df <- data.frame()
  for (i in 1:length(points[,1])) {
    model <- mparse(compartments = compartments, transitions = transitions, u0 = u0, tspan = 1:max(times), gdata = points[i,])
    res <- run(model)
    traj <- trajectory(res)
    trajI <- tidyr::spread(traj[,c('node','time','I')], node, I)
    Iattimes <- trajI[seq_along(trajI[,1])%in%times,]
    aggs <- data.frame(cbind(Iattimes$time, apply(Iattimes[,-1], 1, mean), apply(Iattimes[,-1], 1, function(x) sd(x)/sqrt(nreps)+0.03*(max(x)-min(x))))) %>% setNames(c('time', 'mu', 'sd'))
    if (ev)
      shaped <- c(aggs$mu, aggs$sd) %>% setNames(c(paste0('mu', aggs$time, sep = ""), paste0('EV', aggs$time, sep = "")))
    else
      shaped <- c(aggs$mu) %>% setNames(paste0('mu', aggs$time, sep = ""))
    if (i == 1) {
      out_df <- t(data.frame(shaped))
      rownames(out_df) <- NULL
    }
    else
      out_df <- rbind(out_df, shaped)
  }
  rownames(out_df) <- NULL
  return(out_df)
}

## ----echo = FALSE, message = FALSE--------------------------------------------
library(emulatorr)
library(purrr)
library(SimInf)
library(lhs)

## -----------------------------------------------------------------------------
transitions <- c(
  "S -> beta*I*S/(S+I+E+R) -> E",
  "E -> gamma*E -> I",
  "I -> delta*I -> R",
  "R -> mu*R -> S"
)
compartments <- c("S","E","I","R")
nreps <- 50
u0 <- data.frame(
  S = rep(950, nreps),
  E = rep(0, nreps),
  I = rep(50, nreps),
  R = rep(0, nreps)
)


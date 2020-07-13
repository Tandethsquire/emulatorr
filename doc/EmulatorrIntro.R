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

## -----------------------------------------------------------------------------
params <- c(beta = 0.5, gamma = 0.5, delta = 0.1, mu = 0.1)
model <- mparse(transitions = transitions, compartments = compartments,
                u0 = u0, gdata = params, tspan = 1:60)
result = run(model)
plot(result, compartments = c("E","I","R"))

## -----------------------------------------------------------------------------
points <- expand.grid(list(beta = c(0.4, 0.6),
                           gamma = c(0.4, 0.6),
                           delta = c(0.05, 0.15),
                           mu = c(0.05, 0.15)
))
outputs <- data.frame(cbind(points, getOutputs(points, seq(10,30,by=5))))
head(outputs)

## -----------------------------------------------------------------------------
ranges <- list(
  beta = c(0.2, 0.8),
  gamma = c(0.2, 1),
  delta = c(0.1, 0.5),
  mu = c(0.1, 0.5)
)
pts <- 2*(maximinLHS(80, 4)-1/2)
r_centers <- map_dbl(ranges, ~(.[2]+.[1])/2)
r_scales <- map_dbl(ranges, ~(.[2]-.[1])/2)
pts <- data.frame(t(apply(pts, 1, function(x) x*r_scales + r_centers)))
head(pts)

## -----------------------------------------------------------------------------
wave0 <- data.frame(cbind(pts,getOutputs(pts, seq(10,30,by=5)))) %>% setNames(c(names(ranges), paste0("I",seq(10,30,by=5)),paste0('ev',seq(10,30,by=5))))
head(wave0)

## -----------------------------------------------------------------------------
samples <- sample(nrow(wave0), 40)
train0 <- wave0[samples,1:9]
valid0 <- wave0[!seq_along(wave0[,1])%in%samples,1:9]
output_names <- paste0("I",seq(10,30,by=5))
ems0 <- emulator_from_data(train0, output_names, ranges, quadratic = TRUE)

## -----------------------------------------------------------------------------
delts <- apply(wave0[10:ncol(wave0)], 2, mean)/map_dbl(ems0, ~.$u_sigma)
ems0 <- emulator_from_data(train0, output_names, ranges, quadratic = TRUE, deltas = delts)
ems0[[1]]


## -----------------------------------------------------------------------------
for (i in 1:length(ems0)) ems0[[i]]$output_name <- output_names[i]
names(ems0) <- output_names
emulator_plot(ems0)
emulator_plot(ems0, 'sd')

## -----------------------------------------------------------------------------
wave1 <- map(seq_along(ems0), ~ems0[[.]]$adjust(train0, output_names[[.]]))

## -----------------------------------------------------------------------------
names(wave1) <- output_names
emulator_plot(wave1)
emulator_plot(wave1, var = 'sd')

## -----------------------------------------------------------------------------
em_evals <- wave1$I10$get_exp(train0[,names(ranges)])
all(abs(em_evals - train0$I10) < 10^(-12))
all(wave1$I10$get_cov(train0[,names(ranges)]) < 10^(-12))

## -----------------------------------------------------------------------------
targets = list(
  I10 = list(val = 240, sigma = 25.27),
  I15 = list(val = 396, sigma = 40.99),
  I20 = list(val = 453, sigma = 46.48),
  I25 = list(val = 428, sigma = 43.98),
  I30 = list(val = 392, sigma = 40.30)
)

## -----------------------------------------------------------------------------
emulator_plot(wave1[[1]], 'imp', targets = targets[[1]])

## -----------------------------------------------------------------------------
which_invalid <- validation_diagnostics(wave1, valid0, output_names, targets = targets)
which_invalid

## -----------------------------------------------------------------------------
vp <- validation_pairs(wave1, valid0, targets)

## -----------------------------------------------------------------------------
vp2 <- validation_pairs(wave1, valid0, targets, n=length(wave1))

## -----------------------------------------------------------------------------
sp1 <- space_removed(wave1, valid0, targets)
sp2 <- space_removed(wave1, valid0, targets, u_mod = c(0.75, 1, 1.25), modified = 'var', n_points = 5)
sp3 <- space_removed(wave1, valid0, targets, modified = 'corr', n_points = 5)

## -----------------------------------------------------------------------------
map_dbl(ems0, ~summary(.$model)$adj.r.squared)

## -----------------------------------------------------------------------------
new_points <- generate_new_runs(wave1, ranges, n_points = 120, z = targets)
plot(new_points, pch = 16, cex = 0.5)

## -----------------------------------------------------------------------------
wave_points(list(wave0, new_points), in_names = names(ranges))

## -----------------------------------------------------------------------------
next_wave <- getOutputs(new_points, seq(10,30,by=5))

## -----------------------------------------------------------------------------
w1points <- data.frame(cbind(new_points,next_wave)) %>% setNames(c(names(ranges), paste0("I",seq(10,30,by=5)), paste0("EV",seq(10,30,by=5))))
all_points <- list(wave0[1:9], w1points[1:9])
simulator_plot(all_points, targets)

## -----------------------------------------------------------------------------
sampling <- sample(nrow(w1points), 40)
train1 <- w1points[sampling,1:9]
valid1 <- w1points[!seq_along(w1points[,1])%in%sampling,1:9]
new_ranges <- map(names(ranges), ~c(min(w1points[,.]), max(w1points[,.]))) %>% setNames(names(ranges))
ems1 <- emulator_from_data(train1, output_names, new_ranges, quadratic = T)
deltas <- apply(w1points[,10:14], 2, mean)/map_dbl(ems1, ~.$u_sigma)
ems1 <- emulator_from_data(train1, output_names, new_ranges, deltas = deltas, quadratic = TRUE)
for (i in 1:length(ems1)) ems1[[i]]$output_name <- output_names[i]
wave2 <- map(seq_along(ems1), ~ems1[[.]]$adjust(train1, output_names[[.]]))
names(wave2) <- output_names

## -----------------------------------------------------------------------------
all_waves <- c(wave1, wave2)
all_targets <- c(targets, targets)
emulator_plot(all_waves, var = 'maximp', targets = all_targets)

## -----------------------------------------------------------------------------
test_full_wave <- full_wave(train0, valid0, ranges, output_names, targets, 120, sample_method = 'importance')
wave_points(list(new_points, test_full_wave$next_sample), names(ranges))


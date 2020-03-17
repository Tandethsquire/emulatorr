#' Sample SIR data
#'
#' A small dataset containing points generated using the Gillespie algorithm.
#' The SIR model contains three input parameters, and generates three output
#' parameters. The initial populations are 950 susceptible (S), 50 infected (I),
#' and 0 recovered (R). The final values are taken at time t=20.
#'
#' @format A data frame with 30 rows and 6 variables:
#' \describe{
#'   \item{aSI}{Infection: transition rate from S to I}
#'   \item{aIR}{Recovery: transition rate from I to R}
#'   \item{aSR}{Immunisation: transition rate from S to R}
#'   \item{nS}{Final number of S}
#'   \item{nI}{Final number of I}
#'   \item{nR}{Final number of R}
#' }
"GillespieSIR"

#' Sample SIR validation data
#'
#' A small dataset containing points generated using the Gillespie algorithm.
#' Very similar to \code{\link{GillespieSIR}}, slightly larger in size.
#'
#' @format A data frame with 60 rows and 6 variables:
#' \describe{
#'   \item{aSI}{Infection: transition rate from S to I}
#'   \item{aIR}{Recovery: transition rate from I to R}
#'   \item{aSR}{Immunisation: transition rate from S to R}
#'   \item{nS}{Final number of S}
#'   \item{nI}{Final number of I}
#'   \item{nR}{Final number of R}
#' }
"GillespieValidation"

#' Sample Implausibility Data
#'
#' A dataset containing 1000 points from the region bounded by
#' [0.1, 0.8], [0, 0.5], [0, 0.05] for aSI, aIR and aSR respectively.
#' Implausibility has been calculated (for emulators trained on the
#' \code{\link{GillespieSIR}} dataset) for each of the outputs nS, nI, nR,
#' and the 2nd maximum implausibility is included.
#' The target values used in calculating implausibility were:
#' nS: 281 (sigma 10.43);
#' nI: 30 (sigma 11.16);
#' nR: 689 (sigma 14.32)
#'
#' @format A data frame with 1000 rows and 7 variables:
#' \describe{
#'   \item{aSI}{Infection: transition rate from S to I}
#'   \item{aIR}{Recovery: transition rate from I to R}
#'   \item{aSR}{Immunisation: transition rate from S to R}
#'   \item{nS}{Implausibility for nS}
#'   \item{nI}{Implausibility for nI}
#'   \item{nR}{Implausibility for nR}
#'   \item{I2}{Second-maximum implausibility}
#' }
"GillespieImplausibility"

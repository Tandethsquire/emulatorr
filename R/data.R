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
#' and the maximum implausibility is included.
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
#'   \item{I}{Maximum implausibility}
#' }
"GillespieImplausibility"

#' Sample Multi-wave Emulators
#'
#' An rda object containing three waves of emulation on the
#' Gillespie SIR model.
#'
#' @format A list containing \code{\link{Emulator}} objects:
#' \describe{
#'   \item{Wave 1}{Emulators trained on GillespieSIR to generate wave 2 points}
#'   \item{Wave 2}{Emulators trained on the results of the above wave 2 points}
#'   \item{Wave 3}{Emulators trained on the results of the wave 3 points}
#' }
"GillespieMultiWaveEmulators"

#' Sample Multi-wave Results
#'
#' An rda object containing four data.frames: an initial set of points
#' given by \code{GillespieSIR} and \code{GillespieValidation}, and
#' the 90 points generated at each of three subsequent waves. The trained
#' emulators are provided in \code{\link{GillespieMultiWaveEmulators}}.
#'
#' @format A list of data.frame objects:
#' \describe{
#'   \item{Wave 0}{The initial points used in other examples}
#'   \item{Wave 1}{Points generated from the wave 1 emulators}
#'   \item{Wave 2}{Points generated from the wave 2 emulators}
#'   \item{Wave 3}{Points generated from the wave 3 emulators}
#' }
"GillespieMultiWaveData"

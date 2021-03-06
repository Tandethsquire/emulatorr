% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostics.R
\name{validation_pairs}
\alias{validation_pairs}
\title{Validation Set Comparisons and Implausibility}
\usage{
validation_pairs(ems, validation_points, z, orig_ranges, cb = FALSE, ...)
}
\arguments{
\item{ems}{The list of trained emulators}

\item{validation_points}{The validation set to be plotted}

\item{z}{The target values for each emulated output}

\item{orig_ranges}{The original ranges for the input parameters (if desired)}

\item{cb}{Should a colourblind-friendly palette be used for plots? Default: FALSE}

\item{...}{Any additional parameters to be passed to internal functions.}
}
\value{
A data.frame containing the validation points, with goodness-of-fit and implausibility.
}
\description{
Creates pairs plots on the set of validation points.
}
\details{
Plots are organised as:

a) Emulated vs Simulator Output (lower diagonal). The emulator outputs are compared
against the simulator outputs. Points whose emulated output lies outside the 3-sigma
region of the simulated output are coloured red; those inside are coloured green; a
gradient between the two extremes indicates goodness-of-fit;

b) Implausibility (upper diagonal). The implausibility for each point is calculated,
using the same colour scaling as the lower diagonal.
}
\examples{
ems <- emulator_from_data(GillespieSIR, c('nS','nI','nR'),
 ranges = list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05)),
 quadratic = TRUE)
targets <- list(
 list(val = 281, sigma = 10.43),
 list(val = 30, sigma = 11.16),
 list(val = 689, sigma = 14.32)
)
validation_pairs(ems, GillespieValidation, targets)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derivative_functions.R
\name{directional_fit}
\alias{directional_fit}
\title{Derivative inner product (EXPERIMENTAL)}
\usage{
directional_fit(em, x, v, ...)
}
\arguments{
\item{em}{The emulator in question}

\item{x}{The point in input space from which we want to consider the derivative}

\item{v}{Any direction in the d-dimensional space considered by the emulator}

\item{...}{Additional arguments to be passed along (e.g. local.var to the emulator)}
}
\value{
A 2-vector, consisting of the lower and upper (3-sigma) bounds for the inner product.
}
\description{
Find the (uncertainty-modified) inner product between the derivative at a point \code{x} and
a proposal direction \code{v}.
}
\details{
Given a point \code{x} and a direction \code{v}, we find the overlap between E[f'(x)]
and \code{v}. The emulated derivative has uncertainty
associated with it: this variance is given by v %*% Var[f'(x)] %*% v. Both of the
emulator quantities are obtained from \code{get_deriv_info}.

This function is concerned with ascertaining whether a direction is oriented in the direction
of the emulator gradient. It allows for a consideration of 'emulated gradient descent'.
}

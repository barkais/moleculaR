% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_vibrations.R
\name{ring.vib.df.multi}
\alias{ring.vib.df.multi}
\title{Get a ring's characteristic vibrations frequencies}
\usage{
ring.vib.df.multi(inputs_vector)
}
\arguments{
\item{inputs_vector}{vector of characters of 6 atoms, ordered as follows:
para (opposite of primary - doesn't matter, but must be consistent),
two ortho atoms to primary and two meta atoms to primary.
 See full manual for a detailed explanation.}
}
\value{
A data frame with stretching frequencies
}
\description{
used on a set of molecules, applicable for as many rings wanted
}

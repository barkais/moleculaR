% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_cube_steRimol.R
\name{cube.steRimol.df}
\alias{cube.steRimol.df}
\title{Compute cube based sterimol parameters}
\usage{
cube.steRimol.df(
  coor.atoms,
  only.sub = T,
  drop = NULL,
  plot = F,
  degrees = 90,
  isovalue = 0.003
)
}
\arguments{
\item{coor.atoms}{primary axis, a two atom character}

\item{only.sub}{if TRUE (default) will account only for atoms directly bonded}

\item{drop}{numeric value of an atom, from which onward atoms will be dropped from calculation}

\item{plot}{if TRUE (default) will generate a 3D plot of molecule, with the
density surface and sterimol arrows.}

\item{degrees}{define the rotational scan, defaults to 90 (same as 360)}

\item{isovalue}{what minimal density isovalue should be considered, defaults to 0.003}
}
\value{
sterimol values for a single molecule
}
\description{
The function expects to be executed in a working directory
with cube files. It is recommended
}

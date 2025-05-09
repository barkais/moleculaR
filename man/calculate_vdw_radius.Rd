% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_wVDW_steRimol.R
\name{calculate_vdw_radius}
\alias{calculate_vdw_radius}
\title{Calculate Van der Waals Radius from Non-bonded Atomic Distances}
\usage{
calculate_vdw_radius(
  xyz_data,
  atom_index,
  min_distance = 1,
  max_distance = 5,
  bonds
)
}
\arguments{
\item{xyz_data}{Data frame containing atomic coordinates with columns: index, x, y, z}

\item{atom_index}{Integer specifying the index of the atom to analyze}

\item{min_distance}{Numeric minimum distance threshold in Angstroms (default: 1.0)}

\item{max_distance}{Numeric maximum distance threshold in Angstroms (default: 5.0)}

\item{bonds}{Data frame containing bonding information with columns V1 and V2}
}
\value{
Numeric value of the estimated van der Waals radius in Angstroms
}
\description{
Estimates the van der Waals radius of an atom by analyzing its closest non-bonded contact.
The method assumes that the shortest non-bonded distance represents approximately twice 
the van der Waals radius.
}
\details{
The function excludes bonded atoms and finds the closest non-bonded contact within 
the specified distance range. The van der Waals radius is calculated as half the 
shortest non-bonded distance.
}

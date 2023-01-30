####### ----------------------------------------------------#####
####### -----------------Utility Functions------------------#####
####### ----------------------------------------------------#####

#' Compute angles and dihedrals
#'
#' This function acts on a read feather file
#' Input is a character of either 3 or 4 atoms with
#' atoms separated by a space (e.g. '1 2 3' or '1 2 3 4')
#'
#' @param atoms a character
#' @keywords internal
#' @return 1x1 data frame of the angle/dihedral value in degrees
mol.angles <- function(atoms) {
  atoms.vec <- strsplit(atoms, " ")
  unlisted.atoms <- unlist(atoms.vec)
  numeric.atoms <- as.numeric(unlisted.atoms)
  ang.or.dih <- length(numeric.atoms)
  if (length(numeric.atoms) == 3) {
    numeric.atoms <- c(numeric.atoms[1],
                       rep(numeric.atoms[2], 2),
                       numeric.atoms[3])
  }
  molecule <- list.files(pattern = '.xyz')
  angle.df <- data.frame(matrix(ncol = 1, nrow = 1))
  colnames(angle.df)[1] <- ifelse(ang.or.dih == 3,
    paste('Angle(', stringr::str_replace_all(atoms, ' ', '_'),')', sep = ''),
    paste('Dihedral(', stringr::str_replace_all(atoms, ' ', '_'),')', sep = ''))
  xyz <- data.table::fread(molecule)
  xyz <- xyz[, -1]
  if (ang.or.dih == 3) {
    first.bond <- xyz[numeric.atoms[1], ] - xyz[numeric.atoms[2], ]
    second.bond <- xyz[numeric.atoms[4], ] - xyz[numeric.atoms[3], ]
    angle.df[1, 1] <- angle(as.numeric(first.bond),
                            as.numeric(second.bond)) * (180/pi)
  } else {
    first.bond <- xyz[numeric.atoms[1], ] - xyz[numeric.atoms[2], ]
    second.bond <- xyz[numeric.atoms[3], ] - xyz[numeric.atoms[2], ]
    third.bond <-  xyz[numeric.atoms[4], ] - xyz[numeric.atoms[3], ]
    first.cross <- pracma::cross(as.matrix(first.bond),
                                 as.matrix(second.bond))
    second.cross <- pracma::cross(as.matrix(third.bond),
                                  as.matrix(second.bond))
    angle.df[1, 1] <- angle(as.numeric(first.cross),
                            as.numeric(second.cross)) * (180/pi)
  }
  return(angle.df)
}

#' Compute several angles and dihedrals in a single molecule
#'
#' This function acts on a read feather file
#' Input is a character of either 3 or 4 atoms with
#' atoms separated by a space (e.g. '1 2 3' or '1 2 3 4')
#' and sets gathered in a vector (e.g. c('1 2 3', '1 2 3 4'))
#'
#' @param atoms_vector a vector of characters
#' @keywords internal
#' @return 1x(num of sets) data frame of the angle/dihedral values in degrees
mol.angles.set <- function(atoms_vector) {
  ang.list <- list()
  for (i in 1:length(atoms_vector)) {
    ang.list[[match(i, 1:length(atoms_vector))]] <-  mol.angles(atoms_vector[i])
  }
  ang.df <- do.call(cbind, ang.list)
  return(ang.df)
}

#' Compute several angles and dihedrals in a set of molecules
#'
#' This function acts on a read feather file
#' Input is a character of either 3 or 4 atoms with
#' atoms separated by a space (e.g. '1 2 3' or '1 2 3 4')
#' and sets gathered in a vector (e.g. c('1 2 3', '1 2 3 4'))
#'
#' @param atoms_vector a vector of characters
#' @return (num of molecules)x(num of sets) data frame of the angle/dihedral values in degrees
#' @export
mol.angles.df <- function(atoms_vector) {
  molecules <- list.files()
  mol.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    mol.list[[match(molecule, molecules)]] <- mol.angles.set(atoms_vector)
    setwd('..')
  }
  ang.df <- do.call(rbind, mol.list)
  row.names(ang.df) <- molecules
  return(ang.df)
}

#' Compute atom distances in a single molecule
#'
#' This function acts on a read feather file
#' Input is a character of atoms, splited to pairs
#' atoms separated by a space (e.g. '1 2 3 4') are accounted
#' for as the pairs 1-2 and 3-4
#'
#' @param atom_pairs a character
#' @keywords internal
#' @return 1x(num of pairs) data frame of the distances in Å
atoms.distance <- function(atom_pairs) {
  paired <- atom.pairs.num(atom_pairs)
  molecule <- list.files(pattern = '.xyz')
  mag <- function(vector) {
    sqrt(vector[[1]]^2 + vector[[2]]^2 + vector[[3]]^2)
  }
  dist.df <- data.frame(matrix(ncol = length(paired), nrow = 1))
  names(dist.df) <- stringr::str_replace(as.character(paired),'c','Dist')
  xyz <- data.table::fread(molecule)
  dist.list <- data.frame(matrix(ncol = length(paired), nrow = 1))
  for (i in 1:length(paired)) {
    dist.list[i] <- mag(xyz[paired[[i]][1], 2:4] - xyz[paired[[i]][2], 2:4])
  }
  names(dist.list) <- paired
  dist.df[1, ] <- dist.list
  return(dist.df)
}

#' Compute atom distances of a set of molecules
#'
#' This function acts on a read feather file
#' Input is a character of atoms, splited to pairs
#' atoms separated by a space (e.g. '1 2 3 4') are accounted
#' for as the pairs 1-2 and 3-4
#'
#' @param atom_pairs a character
#' @return (num of molecules)x(num of pairs) data frame of the distances in Å
#' @export
atoms.distance.df <- function(atom_pairs) {
  molecules <- list.files()
  mol.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    mol.list[[match(molecule, molecules)]] <- atoms.distance(atom_pairs)
    setwd('..')
  }
  dist.df <- do.call(rbind, mol.list)
  row.names(dist.df) <- molecules
  return(dist.df)
}

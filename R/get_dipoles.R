####### ----------------------------------------------------#####
####### -----------------Utility Functions------------------#####
####### ----------------------------------------------------#####


#   calculates the center of mass of the "basic" structure.
#   used in context of choice to evaluate the dipole moment's components
#   relative to it (com of basic).

#' Compute molecule's center of mass
#' @param coor_atoms 3 atoms character (e.g. '1 2 3')
#' @keywords internal
#' @return center of mass
center.of.mass <- function(coor_atoms) {
  home <- getwd()
  setwd(paste0('../', list.files('../', pattern = "basic")))
  coor.trans(coor_atoms)
  xyz <- data.table::fread(list.files(pattern = "tc.xyz"), header = F,
                           colClasses = c("character",
                                          "numeric",
                                          "numeric",
                                          "numeric"),
                           stringsAsFactors = F)
  suppressMessages(xyz$V1 <- plyr::mapvalues(xyz$V1,
                                             from = atomic_masses$V1,
                                             to = atomic_masses$V2
  ))
  xyz <- data.frame(transform(xyz, V1 = as.numeric(V1)))
  M <- sum(xyz[,1])
  for (i in 1:dim(xyz)[1]) xyz[i,2:4] <- xyz[i,1] * xyz[i,2:4]
  com <- c(sum(xyz[,2]), sum(xyz[,3]), sum(xyz[,4]))
  com <- (1/M)*com
  unlink(list.files(pattern = "_tc"))
  setwd(home)
  return(com)
}

# center.of.substructure - calculates the center of mass of a substructure in the
#   "basic" structure.
#   used in context of choice to evaluate the dipole moment's components
#   relative to it (com of basic).

#' Compute molecule's substructure center of mass
#'
#' @param sub_atoms atoms character (e.g. '4 5 6 7 8 9')
#' @keywords internal
#' @return center of mass of a substructure
center.of.substructure <- function(sub_atoms) {
  atoms <- strsplit(sub_atoms, " ")
  unlisted.atoms <- unlist(atoms)
  numeric.atoms <- as.numeric(unlisted.atoms)
  xyz <- data.table::fread(list.files(pattern = ".xyz"), header = F,
                           colClasses = c("character",
                                          "numeric",
                                          "numeric",
                                          "numeric"),
                           stringsAsFactors = F)
  xyz <- xyz[numeric.atoms, -1]
  com <- c(sum(xyz[,1]), sum(xyz[,2]), sum(xyz[,3]))
  com <- (1/length(numeric.atoms))*com
  return(com)
}

# npa.dipole.subunit - extarcts the npa charge based dipole moment for a
#   substructure of choice.
#   Expects: coor.atoms - character of three atoms, by which to transform
#   the coordinate system (defines origin and y direction);
#   subunit - character of as many atoms wanted.
#   This is a single molecule function, usable but unusual.

#' Compute nbo based dipole moment for a subset of atoms
#' @param coor_atoms  3 atoms character (e.g. '1 2 3')
#' @param subunit   atoms character (e.g. '4 5 6 7 8 9')
#' @keywords internal
#' @return npa dipole moment for a subsunit
npa.dipole.subunit <- function(coor_atoms, subunit) {
  charges <- data.table::fread(list.files(pattern = 'nbo'))
  charges[[1]] <- as.numeric(charges[[1]])
  colnames(charges) <- "npa"
  coor.trans(coor_atoms)
  transformed_coordinates <- data.table::fread(list.files(pattern = "tc.xyz"), header = F,
                                               colClasses = c("character",
                                                              "numeric",
                                                              "numeric",
                                                              "numeric"),
                                               stringsAsFactors = F)
  transformed_coordinates <- transformed_coordinates[, 2:4]
  dip_comp_mat <- data.frame(cbind(transformed_coordinates, charges))
  subatoms <- strsplit(subunit, " ")
  unlisted.subatoms <- unlist(subatoms)
  numeric.subatoms <- as.numeric(unlisted.subatoms)
  dip_comp_mat <- dip_comp_mat[numeric.subatoms,]
  dip_vec <- as.numeric(vector(length = 3))
  for (i in 1:dim(dip_comp_mat)[[1]]) {
    dip_comp_mat[i, 5] <- dip_comp_mat[i, 1] * dip_comp_mat[i, 4]
    dip_comp_mat[i, 6] <- dip_comp_mat[i, 2] * dip_comp_mat[i, 4]
    dip_comp_mat[i, 7] <- dip_comp_mat[i, 3] * dip_comp_mat[i, 4]
    dip_vec[1] <- sum(dip_comp_mat[, 5])
    dip_vec[2] <- sum(dip_comp_mat[, 6])
    dip_vec[3] <- sum(dip_comp_mat[, 7])
  }
  mag <- function(vector) {
    sqrt(vector[[1]]^2 + vector[[2]]^2 + vector[[3]]^2)
  }
  final <- data.frame(dip_vec[1], dip_vec[2], dip_vec[3], mag(dip_vec))
  colnames(final) <- c("dipNPA_x_subunit", "dipNPA_y_subunit", "dipNPA_z_subunit", "totalNPA_subunit")
  unlink(list.files(pattern = '_tc'))
  return(final)
}

# dip.gaussian - pulls and manipulates dipole moment vector from Gaussian log.
#   Expects quite a bit of information (that is also interdependent):
#   coor_atoms - choose a coordinate system, as the default uses the center of
#   of mass of each molecule as the origin. To avoid any changes to the vector,
#   run the command with no arguments.
#   center_of_mass - uses the center of mass of the "basic" structure. Expects
#   coor_atoms at the same time, on order to have a y direction defined.
#   center_of_substructure - uses the center of mass of a common
#   substructure (e.g. a ring that all molecules posses). Expects coor_atoms.
#   sub_atoms- substructure atoms, as a character (e.g. '1 2 3 4 5 6').
#   This is the single molecule function, usable but unusual.

#' Use Gaussain's dipole moment, and perform manipulations to it.
#' @param coor_atoms 3 atoms character (e.g. '1 2 3')
#' @param center_of_mass logical, should use center of mass of the
#' basic structure as origin or not, if TRUE - center_of_substructure
#' must be FALSE
#' @param center_of_substructure logical, should use center of mass of
#' a substructure as origin or not, if TRUE - center_of_mass must be FALSE
#' @param sub_atoms
#' ONLY if center_of_substructure is TRUE - atoms character (e.g. '4 5 6 7 8 9')
#' @keywords internal
#' @return dipole moment
dip.gaussian <- function(coor_atoms = '',
                         center_of_substructure = F,
                         sub_atoms = NULL,
                         center_of_mass = F) {
  dipole <- data.frame(data.table::fread(list.files(pattern = "dipole")))
  mag <- function(vector) {
    sqrt(vector[[1]]^2 + vector[[2]]^2 + vector[[3]]^2)
  }
  names(dipole) <- c('dip_x', 'dip_y', 'dip_z', 'Total')
  if (!(coor_atoms == '')) {
    atoms <- strsplit(coor_atoms, " ")
    unlisted.atoms <- unlist(atoms)
    numeric.atoms <- as.numeric(unlisted.atoms)
    xyz <- data.table::fread(list.files(pattern = ".xyz"), header = F,
                             colClasses = c("character",
                                            "numeric",
                                            "numeric",
                                            "numeric"),
                             stringsAsFactors = F)
    xyz <- xyz[, -1]
    xyz <- sapply(xyz, as.numeric)

    if (length(numeric.atoms) == 4) {
      new_origin <- (xyz[numeric.atoms[[1]], ] + xyz[numeric.atoms[[2]], ])/2
      new_y <- as.numeric((xyz[numeric.atoms[[3]], ] - new_origin) /
                            mag(xyz[numeric.atoms[[3]], ] - new_origin))
      coplane <- as.numeric((xyz[numeric.atoms[[4]], ] - new_origin) /
                              mag(xyz[numeric.atoms[[4]], ] - new_origin))
    } else {
      if (center_of_mass) {
        com <- center.of.mass(coor_atoms)
        new_origin <- com
        new_y <- as.numeric((xyz[numeric.atoms[[2]], ] - new_origin) /
                              mag(xyz[numeric.atoms[[2]], ] - new_origin))
        coplane <- as.numeric((xyz[numeric.atoms[[3]], ] - new_origin) /
                                mag(xyz[numeric.atoms[[3]], ] - new_origin))
      }
      if (center_of_substructure) {
        if (sub_atoms == 'NA') {
          new_origin <- xyz[numeric.atoms[[1]], ]
          new_y <- as.numeric((xyz[numeric.atoms[[2]], ] - new_origin) /
                                mag(xyz[numeric.atoms[[2]], ] - new_origin))
          coplane <- as.numeric((xyz[numeric.atoms[[3]], ] - new_origin) /
                                  mag(xyz[numeric.atoms[[3]], ] - new_origin))
        } else {
          com <- center.of.substructure(sub_atoms)
          new_origin <- com
          new_y <- as.numeric((xyz[numeric.atoms[[2]], ] - new_origin) /
                                mag(xyz[numeric.atoms[[2]], ] - new_origin))
          coplane <- as.numeric((xyz[numeric.atoms[[3]], ] - new_origin) /
                                  mag(xyz[numeric.atoms[[3]], ] - new_origin))
        }
      }
      if (!center_of_substructure & !center_of_mass) {
        new_origin <- xyz[numeric.atoms[[1]], ]
        new_y <- as.numeric((xyz[numeric.atoms[[2]], ] - new_origin) /
                              mag(xyz[numeric.atoms[[2]], ] - new_origin))
        coplane <- as.numeric((xyz[numeric.atoms[[3]], ] - new_origin) /
                                mag(xyz[numeric.atoms[[3]], ] - new_origin))
      }
    }

    cross_y_coplane <- pracma::cross(coplane, new_y)
    coef_mat <- aperm(array(c(
      new_y,
      coplane,
      cross_y_coplane), dim = c(3, 3)))

    angle_new.y_coplane <- angle(coplane, new_y)
    x_ang_new.y <- pi / 2
    cop_ang_x <- angle_new.y_coplane - x_ang_new.y
    result_vec <- c(0, cos(cop_ang_x), 0)
    new_x <- solve(coef_mat, result_vec)
    new_z <- pracma::cross(new_x, new_y)
    new_basis <- aperm(array(c(new_x, new_y, new_z), dim = c(3, 3)))
    dipole[1:3] <- round(aperm(new_basis %*% as.numeric(dipole[1:3])), 4)
  }
  dipole <- lapply(dipole, as.numeric)
  return(dipole)
}

####### ----------------------------------------------------#####
####### -------------------User Functions-------------------#####
####### ----------------------------------------------------#####



# npa.dipole.subunit.df - extarcts the npa charge based dipole moment for a
# substructure of choice.
# Expects: coor.atoms - character of three atoms, by which to transform
# the coordinate system (defines origin and y direction);
# subunit - character of as many atoms wanted.
#' Extract NBO charge based DM of a substructure alone
#' @param coor_atoms 3 or 4 atoms character (e.g. '1 2 3' or '1 2 3 4')
#' @param sub_atoms atoms character (e.g. '4 5 6 7 8 9')
#' @keywords internal
#' @return dataframe with npa based DM of a substructure
npa.dipole.subunit.df <- function(coor_atoms = NA, sub_atoms = NA) {
  molecules <- list.files()
  if (is.na(coor_atoms)) coor_atoms <- readline('Enter atoms - origin atom, y axis atom and xy plane atom: ')
  if (is.na(sub_atoms)) sub_atoms <- readline('Enter subunit atoms: ')
  dipole.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    dipole.list[[match(molecule, molecules)]] <- round(npa.dipole.subunit(coor_atoms, sub_atoms), 5)
    setwd('..')
  }
  dipole.dafr <- data.frame(data.table::rbindlist(dipole.list))
  row.names(dipole.dafr) <- molecules
  return(dipole.dafr)
}

# npa.dipole.subunit.multi - allows executing for several units. This
#   is the preferred function for use, as it allows for the minimal and maximal
#   use cases.
#   Expects a subunits_inputs_vector - a vector of subunit characters (e.g.
#   c('1 2 3 4 5 6', '15 16 17 18 19 20'))
#' Compute NBO based DM for several substructures
#' @param coor_atoms 3 or 4 atoms character (e.g. '1 2 3')
#' @param subunits_inputs_vector vector of atoms characters
#'  (e.g. c('10,11,12,13', '4 5 6 7 8 9'))
#' @return dataframe with npa based DM of a substructure
#' @export
npa.dipole.subunit.multi <- function(coor_atoms, subunits_inputs_vector) {
  multi.df <- lapply(subunits_inputs_vector,
                     function(x) npa.dipole.subunit.df(coor_atoms, x))
  for (i in 1:length(multi.df)) {
    names(multi.df[[i]]) <- paste0(names(multi.df[[i]]), paste0('_', as.character(i)))
  }
  multi.df.result <- data.frame(do.call(cbind, multi.df))
  return(multi.df.result)
}

# dip.gaussian.df
#
#   pulls and manipulates dipole moment vector from Gaussian log.
#   Expects quite a bit of information (that is also interdependent):
#   coor_atoms - choose a coordinate system, as the default uses the center of
#   of mass of each molecule as the origin. To avoid any changes to the vector,
#   run the command with no arguments.
#   center_of_mass - uses the center of mass of the "basic" structure. Expects
#   coor_atoms at the same time, on order to have a y direction defined.
#   center_of_substructure - uses the center of mass of a common
#   substructure (e.g. a ring that all molecules posses). Expects coor_atoms.
#   sub_atoms- substructure atoms, as a character (e.g. '1 2 3 4 5 6').

#' Pulls and manipulates dipole moment vector.
#' @param coor_atoms 3 or 4 atoms character (e.g. '1 2 3')
#' @param center_of_mass logical, should use center of mass of the
#' basic structure as origin or not, if TRUE - center_of_substructure
#' must be FALSE
#' @param center_of_substructure logical, should use center of mass of
#' a substructure as origin or not, if TRUE - center_of_mass must be FALSE
#' @param sub_atoms
#' ONLY if center_of_substructure is TRUE - atoms character (e.g. '4 5 6 7 8 9')
#' @return dataframe with npa based DM
#' @export
dip.gaussian.df <- function(coor_atoms = '',
                            center_of_substructure = F,
                            sub_atoms = NULL,
                            center_of_mass = F) {
  molecules <- list.files()
  dip.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    dip.list[[match(molecule, molecules)]] <- dip.gaussian(coor_atoms,
                                                           center_of_substructure,
                                                           sub_atoms,
                                                           center_of_mass)
    setwd('..')
  }
  dip.dafr <- data.frame(data.table::rbindlist(dip.list, fill = T))
  row.names(dip.dafr) <- molecules
  return(dip.dafr)
}

# Allows executing for several units. This
# is the preferred function for use, as it allows for the minimal and maximal
# use cases.
# Expects a subunits_inputs_vector - a vector of subunit characters (e.g.
# c('1 2 3 4 5 6', '15 16 17 18 19 20'))

#' Pulls and manipulates dipole moment vector.
#' Allows for use of several substructures.
#' @param coor_atoms_vector vector of sets of 3 or 4 atoms character (e.g. c('1 2 3', '4 5 6'))
#' @param center_of_mass logical, should use center of mass of the
#' basic structure as origin or not, if TRUE - center_of_substructure
#' must be FALSE
#' @param center_of_substructure logical, should use center of mass of
#' a substructure as origin or not, if TRUE - center_of_mass must be FALSE
#' @param subunits_inputs_vector
#' ONLY if center_of_substructure is TRUE - vector of
#' atoms character (e.g. c('10 11 12 13', '4 5 6 7 8 9'))
#' @return data frame with multiple or single DMs
#' @export
dip.gaussian.multi <- function(coor_atoms_vector = c(''),
                                  center_of_substructure = F,
                                  subunits_inputs_vector = NULL,
                                  center_of_mass = F) {
  if (!is.null(subunits_inputs_vector) & !center_of_mass)  {
    while (length(subunits_inputs_vector) < length(coor_atoms_vector)) {
      subunits_inputs_vector <- c(subunits_inputs_vector, 'NA')
    }
    multi.df <- lapply(coor_atoms_vector,
                       function(x) {
                         dip.gaussian.df(x,
                                         center_of_substructure,
                                         subunits_inputs_vector[match(x, coor_atoms_vector)],
                                         center_of_mass)
                       }
                      )
  } else {
    multi.df <- lapply(coor_atoms_vector,
                       function(x) {
                         dip.gaussian.df(x,
                                         center_of_substructure,
                                         subunits_inputs_vector,
                                         center_of_mass)
                       }
                      )
  }

  for (i in 1:length(multi.df)) {
    names(multi.df[[i]]) <- paste0(names(multi.df[[i]]), paste0('_', as.character(i)))
  }
  multi.df.result <- data.frame(do.call(cbind, multi.df))
  multi.df.result <- multi.df.result[!duplicated(as.list(multi.df.result))]
  return(multi.df.result)
}



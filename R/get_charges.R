####### ----------------------------------------------------#####
####### -----------------Utility Functions------------------#####
####### ----------------------------------------------------#####

# nbo.info - reads and organizes nbo charges of wanted atoms.
#   Expects a string of numbers, separated by a space (e.g. '1 2 3').
#   Can be used on its own, but it makes more sense to use the .df version

#' Pull NBO charges for specific atoms
#'
#' Choose atoms for which you wish to have NBO charges presented
#' @param atom_index a character of atoms
#' @keywords internal
#' @return A data frame with NBO charges
nbo.info <- function(atom_index) {
  atom_index <- as.numeric(unlist(strsplit(atom_index, " ")))
  nbos.info <- list.files(pattern = "nbo")
  nbos <- data.frame(data.table::fread(nbos.info))
  nbo.mol <- data.table::transpose(data.frame(nbos[atom_index, ]))
  colnames(nbo.mol) <- paste0('NPA_', as.character(atom_index))
  return(nbo.mol)
}

# hirsh.info - reads and organizes Hirshfeld charges of wanted atoms.
#   Expects a string of numbers, separated by a space (e.g. '1 2 3').
#   Can be used on its own, but it makes more sense to use the .df version

#' Pull Hirshfeld charges for specific atoms
#'
#' Choose atoms for which you wish to have Hirshfeld charges presented
#' @param atom_index a character of atoms
#' @keywords internal
#' @return A data frame with Hirshfeld charges
hirsh.info <- function(atom_index) {
  atom_index <- as.numeric(unlist(strsplit(atom_index, " ")))
  hirshs.info <- list.files(pattern = "Hirshfeld.csv")
  hirshs <- data.frame(data.table::fread(hirshs.info))
  hirsh.mol <- data.table::transpose(data.frame(hirshs[atom_index, ]))
  colnames(hirsh.mol) <- paste0('Hirshfeld_', as.character(atom_index))
  return(hirsh.mol)
}

# cm5.info - reads and organizes CM5 charges of wanted atoms.
#   Expects a string of numbers, separated by a space (e.g. '1 2 3').
#   Can be used on its own, but it makes more sense to use the .df version

#' Pull CM5 charges for specific atoms
#'
#' Choose atoms for which you wish to have CM5 charges presented
#' @param atom_index a character of atoms
#' @keywords internal
#' @return A data frame with CM5 charges
cm5.info <- function(atom_index) {
  atom_index <- as.numeric(unlist(strsplit(atom_index, " ")))
  cm5s.info <- list.files(pattern = "CM5.csv")
  cm5s <- data.frame(data.table::fread(cm5s.info))
  cm5.mol <- data.table::transpose(data.frame(cm5s[atom_index, ]))
  colnames(cm5.mol) <- paste0('CM5_', as.character(atom_index))
  return(cm5.mol)
}

####### ----------------------------------------------------#####
####### -------------------User Functions-------------------#####
####### ----------------------------------------------------#####

# nbo.df- reads and organizes nbo charges of wanted atoms, and comoutes
#   differences for requested pairs.
#   If ran without output, will prompt questions. To avoid charge differences
#   answer '-'. Mark difference.indices = '-' to avoid it to beginwith.
#   Expects a string of numbers, separated by a space (e.g. '1 2 3').
#   Can be used on its own, but it makes more sense to use the .df version

#' Pull NBO charges for specific atoms
#'
#' Choose atoms for which you wish to have NBO charges presented.
#' Charge differences are also available.
#' @param atom_indices a character of atoms
#' @param difference_indices a character of atom pairs
#' @return A data frame with NBO charges and differences
#' @export
nbo.df <- function(atom_indices, difference_indices = NA) {
  molecules <- list.dirs(recursive = F, full.names = F)
  nbo.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    nbo.list[[match(molecule, molecules)]] <- nbo.info(atom_indices)
    setwd('..')
  }
  nbo.dafr <- data.table::rbindlist(nbo.list, fill = T)
  if (difference_indices != '-' & !is.na(difference_indices)) {
    pairs.vec <- strsplit(difference_indices, " ")
    unlisted.pvec <- unlist(pairs.vec)
    paired <- split(
      unlisted.pvec,
      ceiling(seq_along(unlisted.pvec) / 2)
    )
    pairs.vec.atoms <- strsplit(atom_indices, " ")
    unlisted.pvec.atoms <- unlist(pairs.vec.atoms)
    clean.names <- stringr::str_remove_all(names(nbo.dafr), 'NPA_')
    nbo.df.diff <- data.frame(matrix(ncol = length(paired), nrow = nrow(nbo.dafr)))
    diff_names <- vector(length = length(paired))
    for (i in 1:length(diff_names)) {
      diff_names[i] <- paste0('diff_NBO_', paired[[i]][1], '_',
                             paired[[i]][2])
    }
    names(nbo.df.diff) <- diff_names
    for (i in 1:ncol(nbo.df.diff)) {
      paired.1 <- nbo.dafr[[which(unique(unlisted.pvec.atoms) == paired[[i]][1])]]
      paired.2 <- nbo.dafr[[which(unique(unlisted.pvec.atoms) == paired[[i]][2])]]
      nbo.df.diff[, i] <- paired.1 - paired.2
    }
    nbo.dafr <- data.frame(cbind(nbo.dafr, nbo.df.diff))
  } else {
    nbo.dafr <- data.frame(nbo.dafr)
  }
  row.names(nbo.dafr) <- molecules
  return(nbo.dafr)
}

# hirsh.df- reads and organizes Hirshfeld charges of wanted atoms, and comoutes
#   differences for requested pairs.
#   If ran without output, will prompt questions. To avoid charge differences
#   answer '-'. Mark difference.indices = '-' to avoid it to begin with.
#   Expects a string of numbers, separated by a space (e.g. '1 2 3').
#   Can be used on its own, but it makes more sense to use the .df version

#' Pull Hirshfeld charges for specific atoms
#'
#' Choose atoms for which you wish to have Hirshfeld charges presented.
#' Charge differences are also available.
#' @param atom_indices a character of atoms
#' @param difference_indices a character of atom pairs
#' @return A data frame with Hirshfeld charges and differences
#' @export
hirsh.df <- function(atom_indices, difference_indices = NA) {
  molecules <- list.dirs(recursive = F, full.names = F)
  hirsh.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    hirsh.list[[match(molecule, molecules)]] <- hirsh.info(atom_indices)
    setwd('..')
  }
  hirsh.dafr <- data.table::rbindlist(hirsh.list, fill = T)
  if (difference_indices != '-' & !is.na(difference_indices)) {
    pairs.vec <- strsplit(difference_indices, " ")
    unlisted.pvec <- unlist(pairs.vec)
    paired <- split(
      unlisted.pvec,
      ceiling(seq_along(unlisted.pvec) / 2)
    )
    pairs.vec.atoms <- strsplit(atom_indices, " ")
    unlisted.pvec.atoms <- unlist(pairs.vec.atoms)
    clean.names <- stringr::str_remove_all(names(hirsh.dafr), 'Hirshfeld_')
    hirsh.df.diff <- data.frame(matrix(ncol = length(paired), nrow = nrow(hirsh.dafr)))
    diff_names <- vector(length = length(paired))
    for (i in 1:length(diff_names)) {
      diff_names[i] <- paste0('diff_Hirshfeld_', paired[[i]][1], '_',
                              paired[[i]][2])
    }
    names(hirsh.df.diff) <- diff_names
    for (i in 1:ncol(hirsh.df.diff)) {
      paired.1 <- hirsh.dafr[[which(unique(unlisted.pvec.atoms) == paired[[i]][1])]]
      paired.2 <- hirsh.dafr[[which(unique(unlisted.pvec.atoms) == paired[[i]][2])]]
      hirsh.df.diff[, i] <- paired.1 - paired.2
    }
    hirsh.dafr <- data.frame(cbind(hirsh.dafr, hirsh.df.diff))
  } else {
    hirsh.dafr <- data.frame(hirsh.dafr)
  }
  row.names(hirsh.dafr) <- molecules
  return(hirsh.dafr)
}

# cm5.df- reads and organizes CM5 charges of wanted atoms, and comoutes
#   differences for requested pairs.
#   If ran without output, will prompt questions. To avoid charge differences
#   answer '-'. Mark difference.indices = '-' to avoid it to begin with.
#   Expects a string of numbers, separated by a space (e.g. '1 2 3').
#   Can be used on its own, but it makes more sense to use the .df version

#' Pull CM5 charges for specific atoms
#'
#' Choose atoms for which you wish to have CM5 charges presented.
#' Charge differences are also available.
#' @param atom_indices a character of atoms
#' @param difference_indices a character of atom pairs
#' @return A data frame with CM5 charges and differences
#' @export
cm5.df <- function(atom_indices, difference_indices = NA) {
  molecules <- list.dirs(recursive = F, full.names = F)
  cm5.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    cm5.list[[match(molecule, molecules)]] <- cm5.info(atom_indices)
    setwd('..')
  }
  cm5.dafr <- data.table::rbindlist(cm5.list, fill = T)
  if (difference_indices != '-' & !is.na(difference_indices)) {
    pairs.vec <- strsplit(difference_indices, " ")
    unlisted.pvec <- unlist(pairs.vec)
    paired <- split(
      unlisted.pvec,
      ceiling(seq_along(unlisted.pvec) / 2)
    )
    pairs.vec.atoms <- strsplit(atom_indices, " ")
    unlisted.pvec.atoms <- unlist(pairs.vec.atoms)
    clean.names <- stringr::str_remove_all(names(cm5.dafr), 'CM5_')
    cm5.df.diff <- data.frame(matrix(ncol = length(paired), nrow = nrow(cm5.dafr)))
    diff_names <- vector(length = length(paired))
    for (i in 1:length(diff_names)) {
      diff_names[i] <- paste0('diff_CM5_', paired[[i]][1], '_',
                              paired[[i]][2])
    }
    names(cm5.df.diff) <- diff_names
    for (i in 1:ncol(cm5.df.diff)) {
      paired.1 <- cm5.dafr[[which(unique(unlisted.pvec.atoms) == paired[[i]][1])]]
      paired.2 <- cm5.dafr[[which(unique(unlisted.pvec.atoms) == paired[[i]][2])]]
      cm5.df.diff[, i] <- paired.1 - paired.2
    }
    cm5.dafr <- data.frame(cbind(cm5.dafr, cm5.df.diff))
  } else {
    cm5.dafr <- data.frame(cm5.dafr)
  }
  row.names(cm5.dafr) <- molecules
  return(cm5.dafr)
}

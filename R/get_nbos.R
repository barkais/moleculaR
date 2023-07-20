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
      diff_names[i] <- paste0('diff_', paired[[i]][1], '_',
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


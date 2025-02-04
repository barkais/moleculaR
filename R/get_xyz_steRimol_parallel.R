####### ----------------------------------------------------#####
####### -----------------Utility Functions------------------#####
####### ----------------------------------------------------#####


#' Verloop's sterimol values from moleculaR's csv files - parallel version (for Mac and Linux)
#'
#' Calculate a single molecule's sterimol values
#'
#' @param coordinates primary axis, a two atom character
#' @param CPK logical, if TRUE will use CPK radii based on Verloop's atom type,
#' if FALSE will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @param plot.B1 should a plot of B1 and B5 box be presented
#' @keywords internal
#' @return sterimol values for a single molecule
steRimol.parallel <- function(coordinates, CPK = T, only_sub = T, drop = NULL, plot.B1 = T) {
  mol <- list.files(pattern = '.xyz')
  steRimol.xyz.parallel(mol, coordinates, CPK, only_sub, drop, plot.B1)
}

#' Verloop's Sterimol Values from XYZ Files - - parallel version (for Mac and Linux)
#'
#' Calculates Sterimol parameters (B1, B5, L) for a given molecule based on atomic coordinates.
#'
#' @param mol A data frame representing the molecular structure from an XYZ file.
#' @param coordinates A character string specifying the primary axis, defined by two atom indices separated by a space.
#' @param CPK Logical; if TRUE (default), uses CPK radii based on Verloop's atom types; if FALSE, uses Pyykkö's covalent radii.
#' @param only_sub Logical; if TRUE (default), only considers atoms directly bonded from the primary axis onward.
#' @param drop Numeric vector; specifies atom indices to exclude from calculations.
#' @param plot.B1 Logical; if TRUE, generates a plot of the B1 and B5 box.
#'
#' @return A data frame containing Sterimol values (B1, B5, L) and related parameters for the given molecule.
#'
#' @import ggplot2
#' @importFrom stringr str_replace str_split
#' @importFrom data.table fread
#' @importFrom tibble rownames_to_column
#' @importFrom plyr mutate
#' @importFrom parallel mclapply detectCores
#'
#' @keywords internal
steRimol.xyz.parallel <- function(mol, coordinates, CPK = T, only_sub = T, drop = NULL, plot.B1 = T) {
  name_of_axis <- stringr::str_replace(coordinates, ' ', '_')
  origin <- as.numeric(unlist(strsplit(coordinates, " "))[[1]])
  direction <- as.numeric(unlist(strsplit(coordinates, " "))[[2]])
  bonds <- extract.connectivity(mol, threshold_distance = 2.12, keep.HB = F)
  names(bonds) <- c('V1', 'V2')
  if (is.na(bonds[bonds$V1 == direction, ][1, 2])) {
    coordinates <- paste(coordinates, as.character(bonds[bonds$V2 == direction, ][1, 1]), sep = " ")
  } else {
    coordinates <- paste(coordinates, as.character(bonds[bonds$V1 == direction, ][1, 2]), sep = " ")
  }
  if (as.numeric(unlist(strsplit(coordinates, " "))[[3]]) == origin) {
    coordinates <- as.numeric(unlist(strsplit(coordinates, " "))[1:2])
    coordinates <- paste(as.character(coordinates[1]), as.character(coordinates[2]), sep = ' ')
    if (is.na(bonds[bonds$V1 == direction, ][1, 2])) {
      coordinates <- paste(coordinates, as.character(bonds[bonds$V2 == direction, ][2, 1]), sep = " ")
    } else {
      coordinates <- paste(coordinates, as.character(bonds[bonds$V1 == direction, ][2, 2]), sep = " ")
    }
  }
  if (unlist(strsplit(coordinates, " "))[3] == 'NA') {
    coordinates <- as.numeric(unlist(strsplit(coordinates, " "))[1:2])
    coordinates <- paste(as.character(coordinates[1]), as.character(coordinates[2]), sep = ' ')
    remove.direction <- bonds[bonds$V2 != direction, ]
    if (!is.na(remove.direction[remove.direction$V1 == origin, ][1, 2])) {
      coordinates <- paste(coordinates, as.character(remove.direction[remove.direction$V1 == origin, ][1, 2]), sep = " ")
    } else {
      coordinates <- paste(coordinates, as.character(remove.direction[remove.direction$V1 == origin, ][1, 2]), sep = " ")
    }
  }
  coor.trans.file(coordinates, mol)
  mag <- function(vector) {
    sqrt(vector[[1]]^2 + vector[[2]]^2)
  }
  substi <- data.table::fread(list.files(pattern = "tc.xyz"))
  unlink(list.files(pattern = "_tc"))
  substi <- tibble::rownames_to_column(substi)
  if (only_sub == T) {
    all_paths <- find_paths_with_nodes(bonds,
                                       as.numeric(unlist(stringr::str_split(coordinates, ' ')))[1],
                                       as.numeric(unlist(stringr::str_split(coordinates, ' ')))[2])
    rlev <- unique(unlist(all_paths))
    if (!is.null(drop)) {
      rlev <- rlev[!rlev %in% drop]
    }
    substi <- substi[substi$rowname %in% rlev, ]
  }
  substi <- plyr::mutate(substi, Radius = rep(0, nrow(substi)))
  if (CPK == T) {
    bonds.covalent <- extract.connectivity(mol, threshold_distance = 2.12, keep.HB = F)
    a.types <- unlist(lapply(1:nrow(substi), function(x) NOB.Atype(x, substi, bonds.covalent)))
    find.radius <- function(atom) cpk.radii$CPK.radius[cpk.radii$A.type == a.types[[atom]]]
    substi$Radius <- unlist(lapply(1:length(a.types), find.radius))
  } else {
    for (i in 1:dim(substi)[1]) {
      substi$Radius[i] <- Radii.Pyykko$V3[Radii.Pyykko$V2 == substi$V1[i]]
    }
  }
  substi <- plyr::mutate(substi, magnitude = mag(substi[, c(3, 5)]))
  substi <- plyr::mutate(substi, Bs = substi$magnitude + substi$Radius)
  substi <- plyr::mutate(substi, L = substi$V3 + substi$Radius)
  
  # Rough scan for B1 and its location along L
  scans <- 90 / 5
  rough.deg.list <- seq(scans, 90, scans)
  rough.scan.results <- do.call(rbind, mclapply(rough.deg.list,
                                                function(row) b1.scanner(row, substi), mc.cores = detectCores() - 1))
  
  # Define region for fine scan of B1 and its location
  back.ang <- rough.deg.list[which(rough.scan.results == min(rough.scan.results))][1] - scans
  front.ang <- rough.deg.list[which(rough.scan.results == min(rough.scan.results))][1] + scans
  fine.deg.list <- seq(back.ang, front.ang, 1)
  fine.scan.results <- do.call(rbind, mclapply(fine.deg.list,
                                               function(row) b1.scanner(row, substi), mc.cores = detectCores() - 1))
  row.names(fine.scan.results) <- fine.deg.list
  b1s <- fine.scan.results
  
  B1 <- min(b1s)
  if (plot.B1 == T) invisible(b1.scanner(as.numeric(row.names(b1s)[which.min(b1s$B1)]), substi, plot.B1 = T))
  B1_B5_angle <- round(b1s$B1_B5_angle[which.min(b1s$B1)])
  B5 <- max(substi$Bs)
  L <- max(substi$L)
  loc.B5 <- min(substi$V3[which(substi$Bs == B5)])
  result <- unique(round(data.frame(B1, B5, L, loc.B5, B1_B5_angle), 2))
  names(result) <- paste0(names(result), '_', name_of_axis)
  return(result)
}

####### ----------------------------------------------------#####
####### -------------------User Functions-------------------#####
####### ----------------------------------------------------#####

#' Verloop's sterimol values from moleculaR's csv files - - parallel version (for Mac and Linux)
#'
#' Calculate a set of molecules' sterimol values
#'
#' @param coordinates primary axis, a two atom character
#' @param CPK logical, if TRUE will use CPK radii based on Verloop's atom type,
#' if FALSE will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @return sterimol values for a set of molecules
#' @export
steRimol.df.parallel <- function(coordinates, CPK = T, only_sub = T, drop = NULL) {
  molecules <- list.dirs(full.names = F, recursive = F)
  steri.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    steri.list[[match(molecule, molecules)]] <- steRimol.parallel(coordinates, CPK, only_sub, drop, F)
    setwd('..')
  }
  steri.dafr <- data.frame(data.table::rbindlist(steri.list))
  row.names(steri.dafr) <- molecules
  return(steri.dafr)
}

#' Verloop's sterimol values from xyz files - parallel version (for Mac and Linux)
#'
#' Calculate a set of molecules' sterimol values
#'
#' @param coordinates primary axis, a two atom character
#' @param CPK logical, if TRUE will use CPK radii based on Verloop's atom type,
#' if FALSE will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @return sterimol values for a set of molecules
#' @export
steRimol.xyz.df.parallel <- function(coordinates, CPK = T, only_sub = T, drop = NULL) {
  molecules <- list.files(pattern = '.xyz')
  steri.list <- list()
  for (molecule in molecules) {
    steri.list[[match(molecule, molecules)]] <- steRimol.xyz.parallel(molecule, coordinates, CPK, only_sub, drop, F)
  }
  steri.dafr <- data.frame(data.table::rbindlist(steri.list))
  row.names(steri.dafr) <- molecules
  return(steri.dafr)
}

# steRimol.multi - - parallel version (for Mac and Linux)
#   Allows the calculation of several steRimol values, using
#   different primary axes (i.e. for different substituents/substructures).
#   Expects coordinates_vector - a vector of primary axes (e.g. c('1 2', '3 4'))

#' Verloop's sterimol values from moleculaR's csv files
#'
#' Calculate multiple sterimol values for a set of molecules
#'
#' @param coordinates_vector primary axes, a vector of two atom characters
#' @param CPK logical, if TRUE will use CPK radii based on Verloop's atom type,
#' if FALSE will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @return sterimol values (multiple)
#' @export
steRimol.multi.parallel <- function(coordinates_vector, CPK = T, only_sub = T, drop = NULL) {
  multi.df <- lapply(coordinates_vector,
                     function(x) {
                       steRimol.df.parallel(x,
                                   CPK,
                                   only_sub,
                                   drop)
                     }
  )
  multi.df.result <- data.frame(do.call(cbind, multi.df))
  multi.df.result <- multi.df.result[!duplicated(as.list(multi.df.result))]
  return(multi.df.result)
}

#' Verloop's sterimol values from xyz files - parallel version (for Mac and Linux)
#'
#' Calculate multiple sterimol values for a set of molecules
#'
#' @param coordinates_vector primary axes, a vector of two atom characters
#' @param CPK logical, if TRUE will use CPK radii based on Verloop's atom type,
#' if FALSE will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @return sterimol values (multiple)
#' @export
steRimol.xyz.multi.parallel <- function(coordinates_vector, CPK = T, only_sub = T, drop = NULL) {
  multi.df <- lapply(coordinates_vector,
                     function(x) {
                       steRimol.xyz.df.parallel(x,
                                       CPK,
                                       only_sub,
                                       drop)
                     }
  )
  # for (i in 1:length(multi.df)) {
  #   names(multi.df[[i]]) <- paste0(names(multi.df[[i]]), paste0('_', as.character(i)))
  # }
  multi.df.result <- data.frame(do.call(cbind, multi.df))
  multi.df.result <- multi.df.result[!duplicated(as.list(multi.df.result))]
  return(multi.df.result)
}
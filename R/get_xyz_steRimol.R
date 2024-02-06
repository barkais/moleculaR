####### ----------------------------------------------------#####
####### -----------------Utility Functions------------------#####
####### ----------------------------------------------------#####

# b1.scanner - scans for the minimum distance to the rectangular tangent, at
#   a given angle i.

#' scans for B1 vectors
#'
#' @param i degree
#' @param substi truncated xyz
#' @keywords internal
#' @return list of B1 vectors
b1.scanner <- function(i, substi) {
  plane <- substi[, c(3, 5)]
  rot.mat <- matrix(ncol = 2, nrow = 2)
  co.deg <- cos(i * (pi / 180))
  si.deg <- sin(i * (pi / 180))
  rot.mat[1, ] <- c(co.deg, -1 * si.deg)
  rot.mat[2, ] <- c(si.deg, co.deg)
  ind <- NULL
  tc <- matrix(nrow = dim(plane)[[1]], ncol = 2)
  for (i in 1:dim(plane)[[1]]) {
    tc[i, ] <- aperm(rot.mat %*% as.numeric(plane[i, ]))
  }
  tc <- round(tc, 3)
  avs <- abs(c(
    max(tc[, 1]),
    min(tc[, 1]),
    max(tc[, 2]),
    min(tc[, 2])
  ))
  if (min(avs) == 0) {
    if (any(which(avs == min(avs)) %in% c(1, 2))) {
      tc <- round(tc, 1)
      B1 <- max(substi$Radius[which(tc[, 1] == 0, arr.ind = T)])
      B1.loc <- substi$L[which.max(substi$Radius[which(tc[, 1] == 0, arr.ind = T)])]
      return(c(B1, B1.loc))
    }
    if (any(which(avs == min(avs)) %in% c(3, 4))) {
      tc <- round(tc, 1)
      B1 <- max(substi$Radius[which(tc[, 2] == 0, arr.ind = T)])
      B1.loc <- substi$L[which.max(substi$Radius[which(tc[, 2] == 0, arr.ind = T)])]
      return(c(B1, B1.loc))
    }
  } else {
    if (any(which(avs == min(avs)) %in% c(1, 2))) {
      ind <- (which(abs(tc[, 1]) == min(avs), arr.ind = T))[1]
      against <- vector()
      against.loc <- vector()
      if (any(tc[ind, 1] < 0)) {
        new.ind <- which(tc[, 1] == min(tc[, 1]), arr.ind = T)
        tc <- plyr::mutate(as.data.frame(tc), "indi" = dplyr::between(tc[, 1], tc[new.ind[1], 1], tc[new.ind[1], 1] + 1))
        tc[,1] <- -tc[,1]
      } else {
        tc <- plyr::mutate(as.data.frame(tc), "indi" = dplyr::between(tc[, 1], tc[ind[1], 1] - 1, tc[ind[1], 1]))
      }
      for (i in 2:dim(tc)[1]) {
        if (tc[i, 3] == T) {
          against <- append(against, tc[i, 1] + substi$Radius[i])
          against.loc <- append(against.loc, substi$L[i])
        }
        B1 <- ifelse(length(against) > 0, max(against), abs(tc[ind[1], 1]) + substi$Radius[ind[1]])
        B1.loc <- ifelse(length(against) > 0, against.loc[which.max(against)], substi$L[ind[1]])
      }
      return(c(B1, B1.loc))
    }
    if (any(which(avs == min(avs)) %in% c(3, 4))) {
      ind <- (which(abs(tc[, 2]) == min(avs), arr.ind = T))[1]
      against <- vector()
      against.loc <- vector()
      if (any(tc[ind, 2] < 0)) {
        new.ind <- which(tc[, 2] == min(tc[, 2]), arr.ind = T)
        tc <- plyr::mutate(as.data.frame(tc), "indi" = dplyr::between(tc[, 2], tc[new.ind[1], 2], tc[new.ind[1], 2] + 1))
        tc[,2] <- -tc[,2]
      } else {
        tc <- plyr::mutate(as.data.frame(tc), "indi" = dplyr::between(tc[, 2], tc[ind[1], 2] - 1, tc[ind[1], 2]))
      }
      for (i in 2:dim(tc)[1]) {
        if (tc[i, 3] == T) {
          against <- append(against, tc[i, 2] + substi$Radius[i])
          against.loc <- append(against.loc, substi$L[i])
        }
        B1 <- ifelse(length(against) > 0, max(against), abs(tc[ind[1], 2]) + substi$Radius[ind[1]])
        B1.loc <- ifelse(length(against) > 0, against.loc[which.max(against)], substi$L[ind[1]])
      }
      return(c(B1, B1.loc))
    }
  }
}

# NOB.Atype - determines the Verloop atom type as a function of number of bonds.
#   Adapted from Paton's version of sterimol.
#   Vital while using CPK = T.
#' determines atom type based on number of bonds
#'
#' @param row row number in substi (atom)
#' @param substi truncated xyz
#' @param bonds bonds data frame
#' @keywords internal
#' @return list of atom types
NOB.Atype <- function(row, substi, bonds) { # number of bonds for an atom of the substituent
  symbol <- substi[row, 2]
  index <- as.numeric(substi[row, 1])
  nob <- nrow(bonds[bonds$`Atom 1` == index | bonds$`Atom 2` == index, ])
  if (symbol == 'H') result <- 'H'
  if (symbol == 'P') result <- 'P'
  if (symbol == 'F') result <- 'F'
  if (symbol == 'Cl') result <- 'Cl'
  if (symbol == 'Br') result <- 'Br'
  if (symbol == 'I') result <- 'I'
  if (symbol == 'O') {
    if (nob < 1.5) result <- 'O2'
    if (nob > 1.5) result <- 'O'
  }
  if (symbol == 'S') {
    if (nob < 2.5) result <- 'S'
    if (2.5 < nob & nob < 5.5) result <- 'S4'
    if (nob > 5.5) result <- 'S1'
  }
  if (symbol == 'N') {
    if (nob < 2.5) result <- 'C6/N6'
    if (nob > 2.5) result <- 'N'
  }
  if (symbol == 'C') {
    if (nob < 2.5) result <- 'C3'
    if (nob > 2.5 & nob < 3.5) result <- 'C6/N6'
    if (nob > 3.5) result <- 'C'
  }
  if (!symbol %in% c('H', 'P', 'Cl', 'Br',
                     'I', 'O', 'S', 'N', 'C', 'F')) result <- 'X'
  return(result)
}

# steRimol.xyz - compute sterimol values for a given primary axis (so called
#   coordinates). Allows for the use of two radii models, the default being
#   Pyykko's covalent radii, with CPK as an option.
#   By default, the function uses a graph based serach for atoms that should
#   be excluded (i.e. no a part of the substituent/substructure of interest).
#   This is indicated by the argument only_sub = T. In case this is turned F,
#   the function will use the entire structure. Another option is drop - allows
#   the dropping elements that are a part of the substituent, but users would
#   like to exclude. Uses a graph based "cut" of all atoms that stem from some
#   atom. To use this, user will have to indicate the first atom on the
#   substructure they wish to cut out.
#   Expects coordinates of primary axis as a character (e.g. '1 2').
#   This is a single molecule version. Usable, but rarely used.

#' Verloop's sterimol values from moleculaR's csv files
#'
#' Calculate a single molecule's sterimol values
#'
#' @param coordinates primary axis, a two atom character
#' @param CPK logical, if TRUE will use CPK radii based on Verloop's atom type,
#' if FALSE (default) will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @keywords internal
#' @return sterimol values for a single molecule
steRimol <- function(coordinates, CPK = F, only_sub = T, drop = NULL) {
  mol <- list.files(pattern = '.xyz')
  origin <- as.numeric(unlist(strsplit(coordinates, " "))[[1]])
  direction <- as.numeric(unlist(strsplit(coordinates, " "))[[2]])
  bonds <- extract.connectivity(mol)
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
      coordinates <- paste(coordinates, as.character(remove.direction[remove.direction$V2 == origin, ][1, 1]), sep = " ")
    }
  }
  coor.trans(coordinates)
  mag <- function(vector) {
    sqrt(vector[[1]]^2 + vector[[2]]^2)
  }
  substi <- data.table::fread(list.files(pattern = "tc.xyz"))
  unlink(list.files(pattern = "_tc"))
  substi <- tibble::rownames_to_column(substi)
  if (only_sub == T) {
    G <- igraph::graph.data.frame(bonds, directed = T)
    C <- igraph::all_simple_paths(G, as.character(origin), mode = "all")
    rlev <- vector()
    for (i in 1:length(C)) {
      if (length(C[[i]]) > 2) {
        rlev[[i]] <- stringr::str_extract_all(as.character(C[i]), "`[0-9]{1,3}`")
      }
      if (length(C[[i]]) == 2) {
        two_atoms <- stringr::str_extract_all(as.character(C[i]), "[0-9]{1,3}")[[1]][1:4]
        if (as.character(origin) %in% two_atoms && as.character(direction) %in% two_atoms) {
          rlev[[i]] <- stringr::str_extract_all(as.character(C[i]), "`[0-9]{1,3}`")
        }
      }
    }
    remove.vec <- vector()
    if (!is.null(drop)) {
      for (i in 1:length(rlev)) {
        if (any(as.character(paste('`', drop, '`', sep = '')) %in% rlev[[i]][[1]])) {
          remove.vec <- append(remove.vec, i)
        }
      }
      rlev <- rlev[-remove.vec]
    }
    rlev <- rlev[grep(paste("`", direction, "`", sep = ""), rlev)]
    rlev_atoms <- as.numeric(stringr::str_extract_all(unique(unlist(rlev)), "[0-9]{1,3}"))
    substi <- substi[substi$rowname %in% rlev_atoms, ]
  }
  substi <- plyr::mutate(substi, Radius = rep(0, nrow(substi)))
  if (CPK == T) {
    bonds.covalent <- extract.connectivity(mol, threshold_distance = 1.82)
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
  rough.scan.results <- t(data.frame(lapply(rough.deg.list, function(x) b1.scanner(x, substi))))
  b1s <- rough.scan.results[, 1]
  b1s.loc <- rough.scan.results[, 2]
  
  # Define region for fine scan of B1 and its location
  back.ang <- rough.deg.list[which(b1s == min(b1s))][1] - scans
  front.ang <- rough.deg.list[which(b1s == min(b1s))][1] + scans
  fine.deg.list <- seq(back.ang, front.ang, 1)
  fine.scan.results <- t(data.frame(lapply(fine.deg.list, function(x) b1.scanner(x, substi))))
  b1s <- append(b1s, fine.scan.results[, 1])
  b1s.loc <- append(b1s.loc, fine.scan.results[, 2])
  
  B1 <- min(b1s[b1s >= 0])
  loc.B1 <- max(b1s.loc[which(b1s[b1s >= 0] == min(b1s[b1s >= 0], na.rm = TRUE))])
  B5 <- max(substi$Bs)
  L <- max(substi$L)
  loc.B5 <- min(substi$V3[which(substi$Bs == B5)])
  result <- unique(round(data.frame(B1, B5, L, loc.B5, loc.B1), 2))
  return(result)
}

#' Verloop's sterimol values from xyz files
#'
#' Calculate a single molecule's sterimol values
#'
#' @param coordinates primary axis, a two atom character
#' @param CPK logical, if TRUE will use CPK radii based on Verloop's atom type,
#' if FALSE (default) will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @keywords internal
#' @return sterimol values for a single molecule
steRimol.xyz <- function(mol, coordinates, CPK = F, only_sub = T, drop = NULL) {
  origin <- as.numeric(unlist(strsplit(coordinates, " "))[[1]])
  direction <- as.numeric(unlist(strsplit(coordinates, " "))[[2]])
  bonds <- extract.connectivity(mol)
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
    G <- igraph::graph.data.frame(bonds, directed = T)
    C <- igraph::all_simple_paths(G, as.character(origin), mode = "all")
    rlev <- vector()
    for (i in 1:length(C)) {
      if (length(C[[i]]) > 2) {
        rlev[[i]] <- stringr::str_extract_all(as.character(C[i]), "`[0-9]{1,3}`")
      }
      if (length(C[[i]]) == 2) {
        two_atoms <- stringr::str_extract_all(as.character(C[i]), "[0-9]{1,3}")[[1]][1:4]
        if (as.character(origin) %in% two_atoms && as.character(direction) %in% two_atoms) {
          rlev[[i]] <- stringr::str_extract_all(as.character(C[i]), "`[0-9]{1,3}`")
        }
      }
    }
    remove.vec <- vector()
    if (!is.null(drop)) {
      for (i in 1:length(rlev)) {
        if (any(as.character(paste('`', drop, '`', sep = '')) %in% rlev[[i]][[1]])) {
          remove.vec <- append(remove.vec, i)
        }
      }
      rlev <- rlev[-remove.vec]
    }
    rlev <- rlev[grep(paste("`", direction, "`", sep = ""), rlev)]
    rlev_atoms <- as.numeric(stringr::str_extract_all(unique(unlist(rlev)), "[0-9]{1,3}"))
    substi <- substi[substi$rowname %in% rlev_atoms, ]
  }
  substi <- plyr::mutate(substi, Radius = rep(0, nrow(substi)))
  if (CPK == T) {
    bonds.covalent <- extract.connectivity(mol, threshold_distance = 1.82)
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
  rough.scan.results <- t(data.frame(lapply(rough.deg.list, function(x) b1.scanner(x, substi))))
  b1s <- rough.scan.results[, 1]
  b1s.loc <- rough.scan.results[, 2]
  
  # Define region for fine scan of B1 and its location
  back.ang <- rough.deg.list[which(b1s == min(b1s))][1] - scans
  front.ang <- rough.deg.list[which(b1s == min(b1s))][1] + scans
  fine.deg.list <- seq(back.ang, front.ang, 1)
  fine.scan.results <- t(data.frame(lapply(fine.deg.list, function(x) b1.scanner(x, substi))))
  b1s <- append(b1s, fine.scan.results[, 1])
  b1s.loc <- append(b1s.loc, fine.scan.results[, 2])
  
  B1 <- min(b1s[b1s >= 0])
  loc.B1 <- max(b1s.loc[which(b1s[b1s >= 0] == min(b1s[b1s >= 0], na.rm = TRUE))])
  B5 <- max(substi$Bs)
  L <- max(substi$L)
  loc.B5 <- min(substi$V3[which(substi$Bs == B5)])
  result <- unique(round(data.frame(B1, B5, L, loc.B5, loc.B1), 2))
  return(result)
}

####### ----------------------------------------------------#####
####### -------------------User Functions-------------------#####
####### ----------------------------------------------------#####

# steRimol.xyz.df - compute sterimol values for a given primary axis (so called
#   coordinates). Allows for the use of two radii models, the default being
#   Pyykko's covalent radii, with CPK as an option.
#   By default, the function uses a graph based serach for atoms that should
#   be excluded (i.e. no a part of the substituent/substructure of interest).
#   This is indicated by the argument only_sub = T. In case this is turned F,
#   the function will use the entire structure. Another option is drop - allows
#   the dropping elements that are a part of the substituent, but users would
#   like to exclude. Uses a graph based "cut" of all atoms that stem from some
#   atom. To use this, user will have to indicate the first atom on the
#   substructure they wish to cut out.
#   Expects coordinates of primary axis as a character (e.g. '1 2').

#' Verloop's sterimol values from moleculaR's csv files
#'
#' Calculate a set of molecules' sterimol values
#'
#' @param coordinates primary axis, a two atom character
#' @param CPK logical, if TRUE will use CPK radii based on Verloop's atom type,
#' if FALSE (default) will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @return sterimol values for a set of molecules
#' @export
steRimol.df <- function(coordinates, CPK = F, only_sub = T, drop = NULL) {
  molecules <- list.dirs(full.names = F, recursive = F)
  steri.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    steri.list[[match(molecule, molecules)]] <- steRimol(coordinates, CPK, only_sub, drop)
    setwd('..')
  }
  steri.dafr <- data.frame(data.table::rbindlist(steri.list))
  row.names(steri.dafr) <- molecules
  return(steri.dafr)
}

#' Verloop's sterimol values from xyz files
#'
#' Calculate a set of molecules' sterimol values
#'
#' @param coordinates primary axis, a two atom character
#' @param CPK logical, if TRUE will use CPK radii based on Verloop's atom type,
#' if FALSE (default) will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @return sterimol values for a set of molecules
#' @export
steRimol.xyz.df <- function(coordinates, CPK = F, only_sub = T, drop = NULL) {
  molecules <- list.files(pattern = '.xyz')
  steri.list <- list()
  for (molecule in molecules) {
    steri.list[[match(molecule, molecules)]] <- steRimol.xyz(molecule, coordinates, CPK, only_sub, drop)
  }
  steri.dafr <- data.frame(data.table::rbindlist(steri.list))
  row.names(steri.dafr) <- molecules
  return(steri.dafr)
}

# steRimol.multi - Allows the calculation of several steRimol values, using
#   different primary axes (i.e. for different substituents/substructures).
#   Expects coordinates_vector - a vector of primary axes (e.g. c('1 2', '3 4'))

#' Verloop's sterimol values from moleculaR's csv files
#'
#' Calculate multiple sterimol values for a set of molecules
#'
#' @param coordinates_vector primary axes, a vector of two atom characters
#' @param CPK logical, if TRUE will use CPK radii based on Verloop's atom type,
#' if FALSE (default) will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @return sterimol values (multiple)
#' @export
steRimol.multi <- function(coordinates_vector, CPK = F, only_sub = T, drop = NULL) {
  multi.df <- lapply(coordinates_vector,
                     function(x) {
                       steRimol.df(x,
                                   CPK,
                                   only_sub,
                                   drop)
                     }
  )
  for (i in 1:length(multi.df)) {
    names(multi.df[[i]]) <- paste0(names(multi.df[[i]]), paste0('_', as.character(i)))
  }
  multi.df.result <- data.frame(do.call(cbind, multi.df))
  multi.df.result <- multi.df.result[!duplicated(as.list(multi.df.result))]
  return(multi.df.result)
}

#' Verloop's sterimol values from xyz files
#'
#' Calculate multiple sterimol values for a set of molecules
#'
#' @param coordinates_vector primary axes, a vector of two atom characters
#' @param CPK logical, if TRUE will use CPK radii based on Verloop's atom type,
#' if FALSE (default) will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @return sterimol values (multiple)
#' @export
steRimol.xyz.multi <- function(coordinates_vector, CPK = F, only_sub = T, drop = NULL) {
  multi.df <- lapply(coordinates_vector,
                     function(x) {
                       steRimol.xyz.df(x,
                                       CPK,
                                       only_sub,
                                       drop)
                     }
  )
  for (i in 1:length(multi.df)) {
    names(multi.df[[i]]) <- paste0(names(multi.df[[i]]), paste0('_', as.character(i)))
  }
  multi.df.result <- data.frame(do.call(cbind, multi.df))
  multi.df.result <- multi.df.result[!duplicated(as.list(multi.df.result))]
  return(multi.df.result)
}

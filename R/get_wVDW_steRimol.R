#' Calculate Van der Waals Radius from Non-bonded Atomic Distances
#'
#' @description
#' Estimates the van der Waals radius of an atom by analyzing its closest non-bonded contact.
#' The method assumes that the shortest non-bonded distance represents approximately twice 
#' the van der Waals radius.
#'
#' @param xyz_data Data frame containing atomic coordinates with columns: index, x, y, z
#' @param atom_index Integer specifying the index of the atom to analyze
#' @param min_distance Numeric minimum distance threshold in Angstroms (default: 1.0)
#' @param max_distance Numeric maximum distance threshold in Angstroms (default: 5.0)
#' @param bonds Data frame containing bonding information with columns V1 and V2
#'
#' @return Numeric value of the estimated van der Waals radius in Angstroms
#'
#' @details
#' The function excludes bonded atoms and finds the closest non-bonded contact within 
#' the specified distance range. The van der Waals radius is calculated as half the 
#' shortest non-bonded distance.
#'
#' @export
calculate_vdw_radius <- function(xyz_data,
                                 atom_index,
                                 min_distance = 1.0, 
                                 max_distance = 5.0,
                                 bonds) {
  # Get bonded atoms
  bonded_atoms <- unique(c(
    bonds$V1[bonds$V2 == atom_index],
    bonds$V2[bonds$V1 == atom_index]
  ))
  
  # Find closest non-bonded atom
  closest_dist <- Inf
  closest_atom <- NA
  
  # Function to calculate distance between two 3D points
  calculate_distance <- function(coords1, coords2) {
    sqrt(sum((coords1 - coords2)^2))
  }
  
  for (i in as.numeric(xyz_data$index)) {
    if (i != atom_index && !(i %in% bonded_atoms)) {
      dist <- calculate_distance(
        xyz_data[xyz_data$index == as.character(atom_index), 2:4],
        xyz_data[xyz_data$index == as.character(i), 2:4])
      
      if (dist >= min_distance && dist <= max_distance && dist < closest_dist) {
        closest_dist <- dist
        closest_atom <- i
      }
    }
  }
  
  if (length(dist) == 0) {
    warning("No non-bonded atoms found within the specified distance range")
    return(NA)
  }
  
  # Return half of the minimum non-bonded distance as the vdW radius
  return(dist/2)
}

#' Calculate Charge-Weighted Van der Waals Radius
#'
#' @description
#' Calculates a charge-weighted van der Waals radius using partial atomic charges
#' from available charge models (NBO, Hirshfeld, CM5). The weighting is modulated
#' by a damping factor to prevent extreme variations.
#'
#' @param xyz_data Data frame containing atomic coordinates with columns: index, x, y, z
#' @param atom_index Integer specifying the index of the atom to analyze
#' @param min_distance Numeric minimum distance threshold in Angstroms (default: 1.0)
#' @param max_distance Numeric maximum distance threshold in Angstroms (default: 5.0)
#' @param damping_factor Numeric factor to moderate charge-based weighting (default: 0.7)
#' @param bonds Data frame containing bonding information with columns V1 and V2
#' @param charges List of paths to charge model files (NBO, Hirshfeld, CM5)
#'
#' @return Numeric value of the charge-weighted van der Waals radius in Angstroms
#'
#' @details
#' The function uses available charge models to weight the contribution of each atom
#' to the non-bonded distance. The weighting is centered around 0.5 and dampened to
#' prevent extreme variations. If no charge files are available, falls back to the
#' unweighted calculation.
#'
#' @export
calculate_charge_weighted_vdw <- function(xyz_data,
                                          atom_index,
                                          min_distance = 1.0, 
                                          max_distance = 5.0,
                                          damping_factor = 0.7,
                                          bonds,
                                          charges) {
  # Get bonded atoms
  bonded_atoms <- unique(c(
    bonds$V1[bonds$V2 == atom_index],
    bonds$V2[bonds$V1 == atom_index]
  ))
  
  # Find closest non-bonded atom
  closest_dist <- Inf
  closest_atom <- NA
  
  # Function to calculate distance between two 3D points
  calculate_distance <- function(coords1, coords2) {
    sqrt(sum((coords1 - coords2)^2))
  }
  
  for (i in as.numeric(xyz_data$index)) {
    if (i != atom_index && !(i %in% bonded_atoms)) {
      dist <- calculate_distance(
        xyz_data[xyz_data$index == as.character(atom_index), 2:4],
        xyz_data[xyz_data$index == as.character(i), 2:4])
      
      if (dist >= min_distance && dist <= max_distance && dist < closest_dist) {
        closest_dist <- dist
        closest_atom <- i
      }
    }
  }
  
  if (is.na(closest_atom)) {
    warning("No non-bonded atoms found within the specified distance range")
    return(NA)
  }
  
  # Calculate weighted radii for each charge model
  radii <- list()
  atom_charges <- list()
  
  for (model in names(charges)) {
    q_target <- charges[[model]][atom_index]
    q_closest <- charges[[model]][closest_atom]
    
    # Store charges
    atom_charges[[model]] <- list(
      target = q_target,
      closest = q_closest
    )
    
    # Calculate damped weighting
    raw_weight <- abs(q_closest) / (abs(q_target) + abs(q_closest))
    # Center the weight around 0.5 and dampen its deviation
    damped_weight <- 0.5 + (raw_weight - 0.5) * damping_factor
    
    # Calculate radius
    radii[[model]] <- closest_dist * damped_weight
  }
  
  # Calculate average radius if multiple models are available
  avg_radius <- mean(unlist(radii))
  
  return(#list(
    radius = avg_radius
    # radius_by_model = radii,
    # closest_atom = closest_atom,
    # distance = closest_dist,
    # charges = atom_charges,
    # models_used = names(charge_files)
  )#)
}

#' Calculate Sterimol Parameters Using Derived Van der Waals Radii
#'
#' @description
#' Calculates Sterimol parameters (B1, B5, L) using van der Waals radii derived from
#' molecular geometry and optionally weighted by atomic charges.
#'
#' @param mol Character string specifying path to XYZ file
#' @param coordinates Character string specifying primary axis as "atom1 atom2"
#' @param only_sub Logical; if TRUE (default), only considers atoms bonded to primary axis
#' @param drop Numeric vector of atom indices to exclude from calculation
#' @param plot.B1 Logical; if TRUE, generates visualization of B1 and B5 parameters
#'
#' @return Data frame containing:
#' \describe{
#'   \item{B1}{Minimum width perpendicular to primary axis}
#'   \item{B5}{Maximum width perpendicular to primary axis}
#'   \item{L}{Length along primary axis}
#'   \item{loc.B5}{Position of B5 along primary axis}
#'   \item{B1_B5_angle}{Angle between B1 and B5 vectors}
#' }
#'
#' @details
#' Uses derived van der Waals radii from non-bonded contacts, optionally weighted by
#' atomic charges when available. Supports both single-molecule and batch processing.
#'
#' @examples
#' \dontrun{
#' params <- steRimol.xyz.wVDW("molecule.xyz", "1 2")
#' }
#'
#' @seealso
#' \code{\link{steRimol.wVDW.df}} for batch processing
#' \code{\link{steRimol.wVDW.multi}} for multiple primary axes
#'
#' @export
steRimol.xyz.wVDW <- function(mol, coordinates, only_sub = T, drop = NULL, plot.B1 = T) {
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
  # Prepare data for VDW radii calculation
  xyz_data <- data.frame(row.names(substi), substi[,2:4])
  names(xyz_data) <- c('index','x','y','z')
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
  
  # Look for charge files in current directory
  charge_files <- list()
  if (file.exists("nbo.csv")) charge_files[["nbo"]] <- "nbo.csv"
  if (file.exists("Hirshfeld.csv")) charge_files[["Hirshfeld"]] <- "Hirshfeld.csv"
  if (file.exists("CM5.csv")) charge_files[["CM5"]] <- "CM5.csv"
  
  # Check if any charge files were found
  if (length(charge_files) == 0) {
    Radii <- unlist(lapply(as.numeric(unlist(substi[, 1])),
                           function(x) calculate_vdw_radius(xyz_data,
                                                             x, 
                                                             bonds = bonds)))
  } else {
    # Read charge files
    charges <- list()
    for (model in names(charge_files)) {
      charges[[model]] <- data.table::fread(charge_files[[model]])[[1]]
    }
    Radii <- unlist(lapply(as.numeric(unlist(substi[, 1])),
                           function(x) calculate_charge_weighted_vdw(xyz_data,
                                                                     x,
                                                                     damping_factor = 0.7, 
                                                                     bonds = bonds, 
                                                                     charges = charges)))
  }
  
  
  substi <- data.frame(substi, Radii)
  names(substi)[ncol(substi)] <- 'Radius'
  substi <- plyr::mutate(substi, magnitude = mag(substi[, c(3, 5)]))
  substi <- plyr::mutate(substi, Bs = substi$magnitude + substi$Radius)
  substi <- plyr::mutate(substi, L = substi$V3 + substi$Radius)
  
  # Rough scan for B1 and its location along L
  scans <- 90 / 5
  rough.deg.list <- seq(scans, 90, scans)
  rough.scan.results <- do.call(rbind, lapply(rough.deg.list, function(row) b1.scanner(row, substi)))
  
  # Define region for fine scan of B1 and its location
  back.ang <- rough.deg.list[which(rough.scan.results == min(rough.scan.results))][1] - scans
  front.ang <- rough.deg.list[which(rough.scan.results == min(rough.scan.results))][1] + scans
  fine.deg.list <- seq(back.ang, front.ang, 1)
  fine.scan.results <- do.call(rbind, lapply(fine.deg.list, function(row) b1.scanner(row, substi)))
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

#' Calculate Sterimol Parameters Using Derived Van der Waals Radii - in moleculaR csv fies
#'
#' @description
#' Calculates Sterimol parameters (B1, B5, L) using van der Waals radii derived from
#' molecular geometry and optionally weighted by atomic charges.
#' @param coordinates Character string specifying primary axis as "atom1 atom2"
#' @param only_sub Logical; if TRUE (default), only considers atoms bonded to primary axis
#' @param drop Numeric vector of atom indices to exclude from calculation
#' @param plot.B1 Logical; if TRUE, generates visualization of B1 and B5 parameters
#'
#' @return Data frame containing:
#' \describe{
#'   \item{B1}{Minimum width perpendicular to primary axis}
#'   \item{B5}{Maximum width perpendicular to primary axis}
#'   \item{L}{Length along primary axis}
#'   \item{loc.B5}{Position of B5 along primary axis}
#'   \item{B1_B5_angle}{Angle between B1 and B5 vectors}
#' }
#'
#' @details
#' Uses derived van der Waals radii from non-bonded contacts, optionally weighted by
#' atomic charges when available. Supports both single-molecule and batch processing.
#'
#' @seealso
#' \code{\link{steRimol.wVDW.df}} for batch processing
#' \code{\link{steRimol.wVDW.multi}} for multiple primary axes
#'
#' @export
steRimol.wVDW <- function(coordinates, only_sub = T, drop = NULL, plot.B1 = T) {
  mol <- list.files(pattern = '.xyz')
  steRimol.xyz.wVDW(mol, coordinates, only_sub, drop, plot.B1)
}

####### ----------------------------------------------------#####
####### -------------------User Functions-------------------#####
####### ----------------------------------------------------#####

#' Batch Calculate Sterimol Parameters with Derived Van der Waals Radii
#'
#' @description
#' Calculates Sterimol parameters for multiple molecules using derived van der Waals radii.
#'
#' @param coordinates Character string specifying primary axis as "atom1 atom2"
#' @param only_sub Logical; if TRUE (default), only considers atoms bonded to primary axis
#' @param drop Numeric vector of atom indices to exclude from calculation
#'
#' @return Data frame containing Sterimol parameters for each molecule
#'
#' @details
#' Processes multiple molecules in batch, using derived van der Waals radii from
#' non-bonded contacts and atomic charges when available.
#'
#' @export
steRimol.wVDW.df <- function(coordinates, only_sub = T, drop = NULL) {
  molecules <- list.dirs(full.names = F, recursive = F)
  steri.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    steri.list[[match(molecule, molecules)]] <- steRimol.wVDW(coordinates, only_sub, drop, F)
    setwd('..')
  }
  steri.dafr <- data.frame(data.table::rbindlist(steri.list))
  row.names(steri.dafr) <- molecules
  return(steri.dafr)
}

#' Batch Calculate Sterimol Parameters with Derived Van der Waals Radii
#'
#' @description
#' Calculates Sterimol parameters for multiple molecules using derived van der Waals radii.
#'
#' @param coordinates Character string specifying primary axis as "atom1 atom2"
#' @param only_sub Logical; if TRUE (default), only considers atoms bonded to primary axis
#' @param drop Numeric vector of atom indices to exclude from calculation
#'
#' @return Data frame containing Sterimol parameters for each molecule
#'
#' @details
#' Processes multiple molecules in batch, using derived van der Waals radii from
#' non-bonded contacts and atomic charges when available.
#'
#' @export
steRimol.xyz.wVDW.df <- function(coordinates, only_sub = T, drop = NULL) {
  molecules <- list.files(pattern = '.xyz')
  steri.list <- list()
  for (molecule in molecules) {
    steri.list[[match(molecule, molecules)]] <- steRimol.xyz.wVDW(molecule, 
                                                                  coordinates, 
                                                                  only_sub, 
                                                                  drop, 
                                                                  F)
  }
  steri.dafr <- data.frame(data.table::rbindlist(steri.list))
  row.names(steri.dafr) <- molecules
  return(steri.dafr)
}

#' Multi-Axis Sterimol Parameter Calculation
#'
#' @description
#' Calculates charge weighted VDW radii Sterimol parameters 
#' for multiple primary axes across multiple molecules.
#'
#' @param coordinates_vector Character vector of primary axes (e.g., c("1 2", "3 4"))
#' @param only_sub Logical; if TRUE (default), only considers atoms bonded to primary axis
#' @param drop Numeric vector of atom indices to exclude from calculation
#'
#' @return Data frame containing Sterimol parameters for each molecule and axis combination
#'
#' @details
#' Enables calculation of Sterimol parameters for multiple substituents or substructures
#' in a set of molecules, using derived van der Waals radii.
#'
#' @export
steRimol.wVDW.multi <- function(coordinates_vector, only_sub = T, drop = NULL) {
  multi.df <- lapply(coordinates_vector,
                     function(x) {
                       steRimol.wVDW.df(x,
                                       only_sub,
                                       drop)
                     }
  )
  multi.df.result <- data.frame(do.call(cbind, multi.df))
  multi.df.result <- multi.df.result[!duplicated(as.list(multi.df.result))]
  return(multi.df.result)
}

#' Multi-Axis Sterimol Parameter Calculation - xyz
#'
#' @description
#' Calculates charge weighted VDW radii Sterimol parameters 
#' for multiple primary axes across multiple molecules.
#'
#' @param coordinates_vector Character vector of primary axes (e.g., c("1 2", "3 4"))
#' @param only_sub Logical; if TRUE (default), only considers atoms bonded to primary axis
#' @param drop Numeric vector of atom indices to exclude from calculation
#'
#' @return Data frame containing Sterimol parameters for each molecule and axis combination
#'
#' @details
#' Enables calculation of Sterimol parameters for multiple substituents or substructures
#' in a set of molecules, using derived van der Waals radii.
#'
#' @export
steRimol.xyz.multi <- function(coordinates_vector, only_sub = T, drop = NULL) {
  multi.df <- lapply(coordinates_vector,
                     function(x) {
                       steRimol.xyz.wVDW.df(x,
                                            only_sub,
                                            drop)
                     }
  )
  multi.df.result <- data.frame(do.call(cbind, multi.df))
  multi.df.result <- multi.df.result[!duplicated(as.list(multi.df.result))]
  return(multi.df.result)
}

####### ----------------------------------------------------#####
####### -------------Global Utility Functions---------------#####
####### ----------------------------------------------------#####

# coor.trans - transforms molecular coordinates systems - runs in molecular folder
#   Expects input of 3 atoms as character, separated by spaces (e.g. '1 2 3')
#   atom 1 - origin, atom 2 - y direction, atom 3 - xy plane. Works directly
#   on xyz files, creates a new file with the suffix _tc added

#' Transform coordinate system
#'
#' Works directly
#' on xyz files, creates a new file with the suffix _tc added
#' @param coor_atoms 3 or 4 atoms as character, separated by spaces (e.g. '1 2 3'),
#' atom 1 (or atoms 1 and 2) - origin, atom 2 (or 3) - y direction, atom 3 (or 4) - xy plane
#' @return a new xyz file with the suffix _tc added
#' @keywords internal
#' @aliases coor.trans
coor.trans <- function(coor_atoms) {
  atoms <- strsplit(coor_atoms, " ")
  unlisted.atoms <- unlist(atoms)
  numeric.atoms <- as.numeric(unlisted.atoms)
  molecule <- list.files(pattern = '.xyz')
  xyz <- data.frame(data.table::fread(molecule, header = F))
  count.0 <- function(x) sum(cumsum(x != 0) == 0)
  if (count.0(colSums(xyz[, 2:4])) >= 2) {
    col.1 <- count.0(xyz[,2])
    col.2 <- count.0(xyz[,3])
    col.3 <- count.0(xyz[,4])
    place.num <- which.max(c(col.1, col.2, col.3))
    new.row <- data.frame(matrix(nrow = 1, ncol = 4))
    new.row[1, 1] <- 'X'
    new.row[1, place.num + 1] <- 1
    new.row[1, -c(1, place.num + 1)] <- 0
    names(new.row) <- names(xyz)
    xyz <- data.frame(rbind(xyz, new.row))
    numeric.atoms[3] <- nrow(xyz)
  }
  mag <- function(vector) {
    sqrt(vector[[1]]^2 + vector[[2]]^2 + vector[[3]]^2)
  }
  if (length(numeric.atoms) == 4) {
    new_origin <- (xyz[numeric.atoms[[1]], 2:4] + xyz[numeric.atoms[[2]], 2:4])/2
    new_y <- as.numeric((xyz[numeric.atoms[[3]], 2:4] - new_origin) /
                          mag(xyz[numeric.atoms[[3]], 2:4] - new_origin))
    coplane <- as.numeric((xyz[numeric.atoms[[4]], 2:4] - new_origin) /
                            mag(xyz[numeric.atoms[[4]], 2:4] - new_origin))
  } else {
    new_origin <- xyz[numeric.atoms[[1]], 2:4]
    new_y <- as.numeric((xyz[numeric.atoms[[2]], 2:4] - new_origin) /
                          mag(xyz[numeric.atoms[[2]], 2:4] - new_origin))
    coplane <- as.numeric((xyz[numeric.atoms[[3]], 2:4] - new_origin) /
                            mag(xyz[numeric.atoms[[3]], 2:4] - new_origin))
  }
  cross_y_coplane <- pracma::cross(coplane, new_y)
  coef_mat <- aperm(array(c(
    new_y,
    coplane,
    cross_y_coplane
  ),
  dim = c(3, 3)
  ))
  angle_new.y_coplane <- angle(coplane, new_y)
  x_ang_new.y <- pi / 2
  cop_ang_x <- angle_new.y_coplane - x_ang_new.y
  result_vec <- c(0, cos(cop_ang_x), 0)
  new_x <- solve(coef_mat, result_vec)
  new_z <- pracma::cross(new_x, new_y)
  new_basis <- aperm(array(c(new_x, new_y, new_z), dim = c(3, 3)))
  new_origin <- xyz[numeric.atoms[[1]], 2:4]
  new_coordinates <- matrix(nrow = dim(xyz)[[1]], ncol = 3)
  transformed_coordinates <- matrix(nrow = dim(xyz)[[1]], ncol = 3)
  for (i in 1:dim(xyz)[[1]]) {
    new_coordinates[i, ] <- as.numeric(xyz[i, 2:4] - new_origin)
    transformed_coordinates[i, ] <- aperm(new_basis %*% new_coordinates[i, ])
  }
  transformed_coordinates <- round(transformed_coordinates, 4)
  elements <- xyz[, 1]
  transformed_coordinates <- cbind(elements, transformed_coordinates)
  if ('X' %in% transformed_coordinates[, 1]) transformed_coordinates <- transformed_coordinates[1:(nrow(xyz) - 1),]
  colnames(transformed_coordinates) <- names(xyz)
  num.atoms <- nrow(transformed_coordinates)
  m <- as.data.frame(matrix(NA, ncol = 4, nrow = 2))
  m[1, 1] <- num.atoms
  m[is.na(m)] <- ""
  names(m) <- names(xyz)
  transformed_coordinates <- rbind(m, transformed_coordinates)
  transformed_coordinates[transformed_coordinates == "0"] <- "0.0"
  new_xyz <- knitr::kable(transformed_coordinates,
                          format = "simple", row.names = F, col.names = NULL
  )
  new_xyz <- new_xyz[-1]
  new_xyz <- new_xyz[-length(new_xyz)]
  write(new_xyz, paste(tools::file_path_sans_ext(molecule), "_tc", ".xyz", sep = ""))
}

#' Transform coordinate system for a file of choice
#'
#' Works directly
#' on xyz files, creates a new file with the suffix _tc added
#' @param coor_atoms 3 or 4 atoms as character, separated by spaces (e.g. '1 2 3' or '1 2 3 4'),
#' atom 1 - origin, atom 2 - y direction, atom 3 - xy plane
#' @param molecule An xyz file to work on
#' @return a new xyz file with the suffix _tc added
#' @aliases coor.trans.file
#' @export
coor.trans.file <- function(coor_atoms, molecule) {
  atoms <- strsplit(coor_atoms, " ")
  unlisted.atoms <- unlist(atoms)
  numeric.atoms <- as.numeric(unlisted.atoms)
  xyz <- data.frame(data.table::fread(molecule, header = F))
  count.0 <- function(x) sum(cumsum(x != 0) == 0)
  if (count.0(colSums(xyz[, 2:4])) >= 2) {
    col.1 <- count.0(xyz[,2])
    col.2 <- count.0(xyz[,3])
    col.3 <- count.0(xyz[,4])
    place.num <- which.max(c(col.1, col.2, col.3))
    new.row <- data.frame(matrix(nrow = 1, ncol = 4))
    new.row[1, 1] <- 'X'
    new.row[1, place.num + 1] <- 1
    new.row[1, -c(1, place.num + 1)] <- 0
    names(new.row) <- names(xyz)
    xyz <- data.frame(rbind(xyz, new.row))
    numeric.atoms[3] <- nrow(xyz)
  }
  mag <- function(vector) {
    sqrt(vector[[1]]^2 + vector[[2]]^2 + vector[[3]]^2)
  }
  if (length(numeric.atoms) == 4) {
    new_origin <- (xyz[numeric.atoms[[1]], 2:4] + xyz[numeric.atoms[[2]], 2:4])/2
    new_y <- as.numeric((xyz[numeric.atoms[[3]], 2:4] - new_origin) /
                          mag(xyz[numeric.atoms[[3]], 2:4] - new_origin))
    coplane <- as.numeric((xyz[numeric.atoms[[4]], 2:4] - new_origin) /
                            mag(xyz[numeric.atoms[[4]], 2:4] - new_origin))
  } else {
    new_origin <- xyz[numeric.atoms[[1]], 2:4]
    new_y <- as.numeric((xyz[numeric.atoms[[2]], 2:4] - new_origin) /
                          mag(xyz[numeric.atoms[[2]], 2:4] - new_origin))
    coplane <- as.numeric((xyz[numeric.atoms[[3]], 2:4] - new_origin) /
                            mag(xyz[numeric.atoms[[3]], 2:4] - new_origin))
  }
  cross_y_coplane <- pracma::cross(coplane, new_y)
  coef_mat <- aperm(array(c(
    new_y,
    coplane,
    cross_y_coplane
  ),
  dim = c(3, 3)
  ))
  angle_new.y_coplane <- angle(coplane, new_y)
  x_ang_new.y <- pi / 2
  cop_ang_x <- angle_new.y_coplane - x_ang_new.y
  result_vec <- c(0, cos(cop_ang_x), 0)
  new_x <- solve(coef_mat, result_vec)
  new_z <- pracma::cross(new_x, new_y)
  new_basis <- aperm(array(c(new_x, new_y, new_z), dim = c(3, 3)))
  new_origin <- xyz[numeric.atoms[[1]], 2:4]
  new_coordinates <- matrix(nrow = dim(xyz)[[1]], ncol = 3)
  transformed_coordinates <- matrix(nrow = dim(xyz)[[1]], ncol = 3)
  for (i in 1:dim(xyz)[[1]]) {
    new_coordinates[i, ] <- as.numeric(xyz[i, 2:4] - new_origin)
    transformed_coordinates[i, ] <- aperm(new_basis %*% new_coordinates[i, ])
  }
  transformed_coordinates <- round(transformed_coordinates, 4)
  elements <- xyz[, 1]
  transformed_coordinates <- cbind(elements, transformed_coordinates)
  if ('X' %in% transformed_coordinates[, 1]) transformed_coordinates <- transformed_coordinates[1:(nrow(xyz) - 1),]
  colnames(transformed_coordinates) <- names(xyz)
  num.atoms <- nrow(transformed_coordinates)
  m <- as.data.frame(matrix(NA, ncol = 4, nrow = 2))
  m[1, 1] <- num.atoms
  m[is.na(m)] <- ""
  names(m) <- names(xyz)
  transformed_coordinates <- rbind(m, transformed_coordinates)
  transformed_coordinates[transformed_coordinates == "0"] <- "0.0"
  new_xyz <- knitr::kable(transformed_coordinates,
                          format = "simple", row.names = F, col.names = NULL
  )
  new_xyz <- new_xyz[-1]
  new_xyz <- new_xyz[-length(new_xyz)]
  write(new_xyz, paste(tools::file_path_sans_ext(molecule), "_tc", ".xyz", sep = ""))
}
# extract.connectivity - runs in molecular folder
#   default bond distance includes H bonds (2.02 A) - change as needed

#' Create a bonding data frame
#'
#' Works directly on xyz files
#' @param xyz_file An xyz file, automatic for cases in which there's only one.
#' If used directly, be sure to input or else a non stable behavior is expected.
#' @param threshold_distance The upper limit of distance counted as bonded.
#'  Default includes H bonding at 2.02 Å.
#' @param keep.HB should keep H bonds or not? (logical)
#' @return A data frame of bonds
#' @aliases extract.connectivity
#' @export
extract.connectivity <- function(xyz_file = list.files(pattern = '.xyz'),
                                 threshold_distance = 2.20, keep.HB = T) {
  # Read in the XYZ file as a data frame
  xyz_data <- read.table(xyz_file, sep = "", header = FALSE, skip = 2)
  xyz_data <- tibble::rowid_to_column(xyz_data)
  # Extract the atomic numbers and coordinates from the data frame
  atomic_numbers <- xyz_data[, 1]
  coordinates <- xyz_data[, 3:5]
  atom_symbols <- xyz_data[, 2]
  # Use the atomic numbers and coordinates to calculate the distances
  # between each pair of atoms
  distances <- dist(coordinates, diag = T)
  df.distances <- reshape2::melt(as.matrix(distances), varnames = c("a1", "a2"))

  # Remove all unnatural non-covalent bonding
  df.distances <- dplyr::mutate(df.distances, first = rep(NA, nrow(df.distances)))
  df.distances <- dplyr::mutate(df.distances, second = rep(NA, nrow(df.distances)))
  for (i in 1:nrow(df.distances)) df.distances$first[i] <- atom_symbols[df.distances$a1[i]]
  for (i in 1:nrow(df.distances)) df.distances$second[i] <- atom_symbols[df.distances$a2[i]]
  remove <- vector()
  for (i in 1:nrow(df.distances)) {
    if (df.distances[i, 4] == 'H' & !df.distances[i, 5] %in% c('N', 'O', 'F')) {
      remove <- append(remove, i)
    }
    if (df.distances[i, 4] == 'H' & df.distances[i, 5] == 'H') remove <- append(remove, i)
  }
  if (length(remove) > 0) df.distances <- df.distances[-remove, ]
  remove <- vector()
  H.num <- c(NULL, NULL)
  for (i in 1:nrow(df.distances)) {
    if ((df.distances[i, 4] == 'H' | df.distances[i, 5] == 'H') &
        (1.5 <= df.distances$value[i] & df.distances$value[i] <= threshold_distance)) {
      if (df.distances[i, 4] == 'H') H.num <- c(df.distances$a1[i], i)
      if (df.distances[i, 5] == 'H') H.num <- c(df.distances$a2[i], i)
      if (!is.null(H.num) & keep.HB == T) {
          examine.df <- df.distances[df.distances$a1 == H.num[1] | df.distances$a2 == H.num[1], ]
          examine.df <- examine.df[examine.df$value <= threshold_distance, ]
          if (any(!examine.df$first %in% c('N', 'O', 'F', 'H') | !examine.df$second %in% c('N', 'O', 'F', 'H'))) {
            remove <- append(remove, H.num[2])
          }
        } 
      }
    if ((df.distances[i, 4] == 'H' | df.distances[i, 5] == 'H') &
        (1.5 <= df.distances$value[i] & df.distances$value[i] <= threshold_distance &
         keep.HB == F)) {
      remove <- append(remove, i)
    }
  }
  if (length(remove) > 0) df.distances <- df.distances[-remove, ]
  row.names(df.distances) <- 1:nrow(df.distances)
  # Extract the atom pairs that are within a certain distance of each other
  # (i.e. the atoms that are bonded to each other)
  bonded_pairs <- which(df.distances$value <= threshold_distance &
                          df.distances$value > 0, arr.ind = T)
  bonded_pairs_df <- unique(as.data.frame(t(apply(df.distances[bonded_pairs,
                                                               1:2],
                                                  1,
                                                  sort))))
  colnames(bonded_pairs_df) <- c("Atom 1", "Atom 2")
  return(bonded_pairs_df)
}

# angle - returns the angle between two vectors (degrees)
#
#' Calculate the angle between two vectors (degrees)
#'
#' @param x the first vector
#' @param y the second vector
#' @keywords internal
#' @return angle in degrees
angle <- function(x,y){
  dot.prod <- x %*% y
  norm.x <- norm(x, type = "2")
  norm.y <- norm(y, type = "2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

#' Change File Names in a Directory
#'
#' @description 
#' Renames files in a directory either by direct file mapping or by pattern replacement.
#'
#' @param dir Character string specifying the directory containing files to rename
#' @param current_names Character vector of current filenames or a pattern to replace
#' @param new_names Character vector of replacement filenames or a replacement pattern
#' @param pattern_mode Logical, if TRUE treats the inputs as patterns to replace rather than exact filenames
#' @param add_extension Logical, whether to automatically append original file extension (only used when pattern_mode=FALSE)
#' 
#' @details
#' The function works in two modes:
#' 1. Direct renaming (pattern_mode=FALSE): Maps specific files to new names
#' 2. Pattern replacement (pattern_mode=TRUE): Replaces occurrences of a pattern in filenames
#'
#' @return Character vector of the new filenames
#'
#' @examples
#' \dontrun{
#' # Direct renaming example
#' name_changer("path/to/dir", 
#'              c("file1.txt", "file2.txt"), 
#'              c("newname1", "newname2"))
#'              
#' # Pattern replacement example
#' name_changer("path/to/dir", 
#'              "CE_No_Water", 
#'              "conformer_3", 
#'              pattern_mode = TRUE)
#' }
#'
#' @export
name_changer <- function(dir, current_names, new_names, pattern_mode = FALSE, add_extension = TRUE) {
  # Store original working directory
  original_dir <- getwd()
  on.exit(setwd(original_dir))
  
  # Change to target directory
  setwd(dir)
  
  # Get existing files
  existing_files <- list.files(recursive = FALSE, full.names = FALSE)
  
  # Pattern replacement mode
  if (pattern_mode) {
    if (length(current_names) != 1 || length(new_names) != 1) {
      stop("In pattern mode, current_names and new_names must be single strings")
    }
    
    # Find files matching the pattern
    pattern_files <- grep(current_names, existing_files, value = TRUE)
    
    if (length(pattern_files) == 0) {
      warning("No files found containing the pattern: ", current_names)
      return(character(0))
    }
    
    # Create new names by replacing the pattern
    new_filenames <- gsub(current_names, new_names, pattern_files)
    
    # Check for name collisions
    conflicts <- new_filenames[new_filenames %in% existing_files & 
                                 new_filenames != pattern_files]
    if (length(conflicts) > 0) {
      warning("Some target filenames already exist and will be skipped: ",
              paste(conflicts, collapse = ", "))
      
      # Remove conflicting renames
      keep <- !(new_filenames %in% conflicts)
      pattern_files <- pattern_files[keep]
      new_filenames <- new_filenames[keep]
    }
    
    # Perform renaming
    for (i in seq_along(pattern_files)) {
      if (!file.exists(new_filenames[i]) || new_filenames[i] == pattern_files[i]) {
        file.rename(pattern_files[i], new_filenames[i])
      }
    }
    
    return(new_filenames)
    
  } else {
    # Direct renaming mode
    # Input validation
    if (length(current_names) != length(new_names)) {
      stop("current_names and new_names must have the same length")
    }
    
    # Extract extensions if needed
    if (add_extension) {
      extensions <- tools::file_ext(current_names)
      new_names <- ifelse(extensions != "",
                          paste0(new_names, ".", extensions),
                          new_names)
    }
    
    # Check for existing files to avoid collisions
    collisions <- new_names[new_names %in% existing_files & 
                              !(new_names %in% current_names)]
    if (length(collisions) > 0) {
      warning("Some target filenames already exist and will be skipped: ",
              paste(collisions, collapse = ", "))
    }
    
    # Perform renaming 
    for (i in seq_along(current_names)) {
      if (current_names[i] %in% existing_files && 
          !file.exists(new_names[i])) {
        file.rename(current_names[i], new_names[i])
      }
    }
    
    return(new_names)
  }
}

# atom.pairs.num - breaks a string of atoms, separated by spaces, into pairs

#' Helper function for manipulating human input
#'
#' Breaks a string of atoms, separated by spaces, into pairs
#' @param atom_pairs A character of atoms (even number of atoms)
#' @keywords internal
#' @return list of numeric atom pairs
atom.pairs.num <- function(atom_pairs) {
  bonds.vec <- strsplit(atom_pairs, " ")
  unlisted.bvec <- unlist(bonds.vec)
  numeric.bvec <- as.numeric(unlisted.bvec)
  paired <- split(
    numeric.bvec,
    ceiling(seq_along(numeric.bvec) / 2)
  )
  return(paired)
}


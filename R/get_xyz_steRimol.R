####### ----------------------------------------------------#####
####### -----------------Utility Functions------------------#####
####### ----------------------------------------------------#####

# b1.scanner - scans for the minimum distance to the rectangular tangent, at
#   a given angle i.

#' scans for B1 vectors
#'
#' @param i degree
#' @param plane a 2D projection of the substituent's xyz, perpendicular to L axis
#' @param plot.B1 should a plot of B1 and B5 box be presented
#' @keywords internal
#' @import ggplot2
#' @importFrom data.table data.table :=
#' @return list of B1 vectors
b1.scanner <- function(i, substi, plot.B1 = F) {
  # Extract the plane coordinates and radii
  plane <- substi[, c(3, 5, 6)]  # x, y, radius columns
  
  # Function to generate points for a circle
  generate_circle_points <- function(center_x, center_y, radius, n_points = 100) {
    theta <- seq(0, 2 * pi, length.out = n_points)
    x <- center_x + radius * cos(theta)
    y <- center_y + radius * sin(theta)
    data.frame(x = x, y = y)
  }
  
  # Generate circle points for each atom
  all_points <- list()
  for(j in 1:nrow(substi)) {
    circle <- generate_circle_points(
      substi$V2[j], 
      substi$V4[j], 
      substi$Radius[j]
    )
    all_points[[j]] <- circle
  }
  
  # Combine all points
  points <- do.call(rbind, all_points)
  
  # Create rotation matrix
  angle <- i * (pi / 180)
  rot.mat <- matrix(
    c(cos(angle), -sin(angle),
      sin(angle), cos(angle)),
    nrow = 2, 
    byrow = TRUE
  )
  
  # Rotate all points
  rotated_points <- as.matrix(points) %*% rot.mat
  # Convert to data frame with proper column names
  rotated_df <- data.frame(x = rotated_points[,1], y = rotated_points[,2])
  
  # Find extremes for B1 calculation
  max_x <- max(rotated_df$x)
  min_x <- min(rotated_df$x)
  max_y <- max(rotated_df$y)
  min_y <- min(rotated_df$y)
  
  # Calculate distances to axes
  distances <- abs(c(max_x, min_x, max_y, min_y))
  min_val <- min(distances)
  min_index <- which.min(distances)
  
  # Get B1 coordinates
  b1_coords <- switch(min_index,
                      c(max_x, 0),    # max_x
                      c(min_x, 0),    # min_x
                      c(0, max_y),    # max_y
                      c(0, min_y)     # min_y
  )
  
  # Find B5 point (maximum distance from origin)
  distances_from_origin <- sqrt(rotated_df$x^2 + rotated_df$y^2)
  b5_index <- which.max(distances_from_origin)
  b5_point <- c(rotated_df$x[b5_index], rotated_df$y[b5_index])
  b5_value <- distances_from_origin[b5_index]
  
  # Calculate angle between B1 and B5
  angle_b1 <- atan2(b1_coords[2], b1_coords[1]) %% (2 * pi)
  angle_b5 <- atan2(b5_point[2], b5_point[1]) %% (2 * pi)
  angle_diff <- abs(angle_b5 - angle_b1)
  if(angle_diff > pi) angle_diff <- 2*pi - angle_diff
  
  # Plot if requested
  if(plot.B1) {
    plot_limits <- max(abs(c(min_x, max_x, min_y, max_y))) * 1.2
    
    # Create plot
    plot <- ggplot2::ggplot() +
      # Add rotated circle points
      ggplot2::geom_point(data = rotated_df, 
                          ggplot2::aes(x = x, y = y),
                          color = "cadetblue2", alpha = 0.5) +
      # Add guide lines
      ggplot2::geom_vline(xintercept = c(min_x, max_x), 
                          color = "darkred", linetype = "dashed") +
      ggplot2::geom_hline(yintercept = c(min_y, max_y), 
                          color = "darkgreen", linetype = "dashed") +
      # Add B1 arrow
      ggplot2::geom_segment(ggplot2::aes(x = 0, y = 0, 
                                         xend = b1_coords[1], 
                                         yend = b1_coords[2]),
                            arrow = ggplot2::arrow(type = "closed", 
                                                   length = ggplot2::unit(0.2, "inches")),
                            color = "#8FBC8F", size = 1.2) +
      # Add B5 arrow
      ggplot2::geom_segment(ggplot2::aes(x = 0, y = 0,
                                         xend = b5_point[1], 
                                         yend = b5_point[2]),
                            arrow = ggplot2::arrow(type = "closed", 
                                                   length = ggplot2::unit(0.2, "inches")),
                            color = "#CD3333", size = 1.2) +
      # Add labels
      ggplot2::geom_label(ggplot2::aes(x = b1_coords[1] * 0.5, 
                                       y = b1_coords[2] * 0.5,
                                       label = sprintf("B1\n%.2f", min_val)),
                          color = "black", fill = "white") +
      ggplot2::geom_label(ggplot2::aes(x = b5_point[1] * 0.5,
                                       y = b5_point[2] * 0.5, 
                                       label = sprintf("B5\n%.2f", b5_value)),
                          color = "black", fill = "white") +
      # Add angle label
      ggplot2::geom_label(ggplot2::aes(x = 0.8 * cos((angle_b1 + angle_b5)/2),
                                       y = 0.8 * sin((angle_b1 + angle_b5)/2),
                                       label = sprintf("%.1f\u00b0", angle_diff * 180/pi)),
                          color = "black", fill = "white") +
      # Customize theme
      ggplot2::theme_minimal() +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::ggtitle(sprintf("Plane rotated %d degrees", i)) +
      ggplot2::xlim(-plot_limits, plot_limits) +
      ggplot2::ylim(-plot_limits, plot_limits)
    
    print(plot)
  }
  
  return(data.frame(B1 = min_val, B1_B5_angle = angle_diff * 180/pi))
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
  if (!symbol %in% c('H', 'P', 'F','Cl', 'Br',
                     'I', 'O', 'S', 'N', 'C')) result <- 'X'
  return(result)
}

# find_paths_with_nodes - find all simple paths that include the 
# primaru steRimol axis. Used in the process of steRimol computation 
#' when only.sub is TRUE.
#' @param bonds bonds dataframe
#' @param node1 primary steRimol atom
#' @param bonds secondary steRimol atom
#' @keywords internal
#' @return all simple paths that include the two atoms
find_paths_with_nodes <- function(bonds, node1, node2) {
  
  # Convert the dataframe to an adjacency matrix
  graph <- matrix(0, nrow = max(max(bonds)), ncol = max(max(bonds)))
  for (i in 1:nrow(bonds)) {
    graph[bonds[i, 1], bonds[i, 2]] <- 1
    graph[bonds[i, 2], bonds[i, 1]] <- 1
  }
  
  all_paths <- list()
  
  # Depth-first search function
  dfs <- function(current_node, current_path) {
    current_path <- c(current_path, current_node)
    
    # If the current path includes both nodes, add it to the list of paths
    if (node1 %in% current_path & node2 %in% current_path) {
      all_paths <<- c(all_paths, list(current_path))
    }
    
    # Recursive DFS for all neighbors
    neighbors <- which(graph[current_node, ] == 1)
    for (neighbor in neighbors) {
      if (!(neighbor %in% current_path)) {
        dfs(neighbor, current_path)
      }
    }
  }
  
  # Start DFS from each node
  for (start_node in node1:node1) {
    dfs(start_node, numeric(0))
  }
  
  return(all_paths)
}

# steRimol - in a moleculaR_csv_files folder -compute sterimol values for a
#   given primary axis (so called coordinates).
#   Allows for the use of two radii models, the default being
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
#' if FALSE will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @param plot.B1 should a plot of B1 and B5 box be presented
#' @keywords internal
#' @return sterimol values for a single molecule
steRimol <- function(coordinates, CPK = T, only_sub = T, drop = NULL, plot.B1 = T) {
  mol <- list.files(pattern = '.xyz')
  steRimol.xyz(mol, coordinates, CPK, only_sub, drop, plot.B1)
}

#' Verloop's Sterimol Values from XYZ Files
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
#' @importFrom dplyr mutate
#'
#' @keywords internal
steRimol.xyz <- function(mol, coordinates, CPK = T, only_sub = T, drop = NULL, plot.B1 = T) {
  name_of_axis <- stringr::str_replace(coordinates, ' ', '_')
  origin <- as.numeric(unlist(strsplit(coordinates, " "))[[1]])
  direction <- as.numeric(unlist(strsplit(coordinates, " "))[[2]])
  bonds <- extract.connectivity(mol, threshold_distance = 2.20, keep.HB = F)
  names(bonds) <- c('V1', 'V2')
  if (is.na(bonds[bonds$V1 == direction, ][1, 2])) {
    coordinates <- paste(coordinates, as.character(bonds[bonds$V2 == direction, ][1, 1]), sep = " ")
  } else {
    coordinates <- paste(coordinates, as.character(bonds[bonds$V1 == direction, ][1, 2]), sep = " ")
  }
  if (unlist(strsplit(coordinates, " "))[[3]] == as.character(origin)) {
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
  substi <- dplyr::mutate(substi, Radius = rep(0, nrow(substi)))
  if (CPK == T) {
    bonds.covalent <- extract.connectivity(mol, threshold_distance = 2.20, keep.HB = F)
    a.types <- unlist(lapply(1:nrow(substi), function(x) NOB.Atype(x, substi, bonds.covalent)))
    find.radius <- function(atom) cpk.radii$CPK.radius[cpk.radii$A.type == a.types[[atom]]]
    substi$Radius <- unlist(lapply(1:length(a.types), find.radius))
  } else {
    for (i in 1:dim(substi)[1]) {
      substi$Radius[i] <- Radii.MMP$V3[Radii.MMP$V2 == substi$V1[i]]
    }
  }
  substi <- dplyr::mutate(substi, magnitude = mag(substi[, c(3, 5)]))
  substi <- dplyr::mutate(substi, Bs = substi$magnitude + substi$Radius)
  substi <- dplyr::mutate(substi, L = substi$V3 + substi$Radius)
  
  # Rough scan for B1 and its location along L
  scans <- 90 / 5
  rough.deg.list <- seq(scans, 90, scans)
  rough.scan.results <- do.call(rbind, lapply(rough.deg.list, function(row) b1.scanner(row, substi)))
  
  # Define region for fine scan of B1 and its location
  back.ang <- rough.deg.list[which(rough.scan.results$B1 == min(rough.scan.results$B1))][1] - scans
  front.ang <- rough.deg.list[which(rough.scan.results$B1 == min(rough.scan.results$B1))][1] + scans
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
#' if FALSE will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @return sterimol values for a set of molecules
#' @export
steRimol.df <- function(coordinates, CPK = T, only_sub = T, drop = NULL) {
  molecules <- list.dirs(full.names = F, recursive = F)
  steri.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    steri.list[[match(molecule, molecules)]] <- steRimol(coordinates, CPK, only_sub, drop, F)
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
#' if FALSE will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @return sterimol values for a set of molecules
#' @export
steRimol.xyz.df <- function(coordinates, CPK = T, only_sub = T, drop = NULL) {
  molecules <- list.files(pattern = '.xyz')
  steri.list <- list()
  for (molecule in molecules) {
    steri.list[[match(molecule, molecules)]] <- steRimol.xyz(molecule, coordinates, CPK, only_sub, drop, F)
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
#' if FALSE will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @return sterimol values (multiple)
#' @export
steRimol.multi <- function(coordinates_vector, CPK = T, only_sub = T, drop = NULL) {
  multi.df <- lapply(coordinates_vector,
                     function(x) {
                       steRimol.df(x,
                                   CPK,
                                   only_sub,
                                   drop)
                     }
  )
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
#' if FALSE will use Pyykko's covalent radii
#' @param only_sub if TRUE (default) will account only for atoms directly bonded
#' from the primary axis onward
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @return sterimol values (multiple)
#' @export
steRimol.xyz.multi <- function(coordinates_vector, CPK = T, only_sub = T, drop = NULL) {
  multi.df <- lapply(coordinates_vector,
                     function(x) {
                       steRimol.xyz.df(x,
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
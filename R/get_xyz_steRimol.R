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
#' @importFrom data.table data.table :=
#' @return list of B1 vectors
b1.scanner <- function(i, substi, plot.B1 = F) {
  plane <- substi[, c(3, 5, 6)]
  generate_circle <- function(center_x, center_y, radius, n_points = 20) {
    theta <- seq(0, 2 * pi, length.out = n_points)
    data.table::data.table(
      x = center_x + radius * cos(theta),
      y = center_y + radius * sin(theta)
    )
  }
  
  # Generate circles for each point
  circle_points <- plane[, generate_circle(V2, V4, Radius), by = .(id = .I)]
  plane <- circle_points[, 2:3]
  rot.mat <- matrix(ncol = 2, nrow = 2)
  co.deg <- cos(i * (pi / 180))
  si.deg <- sin(i * (pi / 180))
  rot.mat[1, ] <- c(co.deg, -1 * si.deg)
  rot.mat[2, ] <- c(si.deg, co.deg)
  
  tc <- matrix(nrow = dim(plane)[[1]], ncol = 2)
  for (row in 1:dim(plane)[[1]]) {
    tc[row, ] <- aperm(rot.mat %*% as.numeric(plane[row, ]))
  }
  tc <- round(tc, 3)
  
  # Calculate distances to lines
  max_x <- max(tc[, 1])
  min_x <- min(tc[, 1])
  max_y <- max(tc[, 2])
  min_y <- min(tc[, 2])
  
  avs <- abs(c(max_x, min_x, max_y, min_y))
  
  min_val <- min(avs)
  min_index <- which.min(avs)
  
  # Determine the coordinates for B1 arrow based on which value was minimum
  b1_coords <- switch(min_index,
                      c(max_x, 0),    # if max_x was minimum
                      c(min_x, 0),    # if min_x was minimum
                      c(0, max_y),    # if max_y was minimum
                      c(0, min_y))    # if min_y was minimum
  
  
  
  
  # B5 for angle between them and for visualization
  mag <- function(vector) {
    sqrt(sum((vector)^2))
  }
  
  b5.point <- tc[which.max(apply(tc, 1, mag)), ]
  b5.value <- mag(b5.point)
  
  
  # Calculate angle between B1 and B5 arrows
  # Using atan2 for correct quadrant handling
  angle_b1 <- atan2(b1_coords[2], b1_coords[1])
  angle_b5 <- atan2(b5.point[2], b5.point[1])    
  
  # Ensure angles are positive
  angle_b1 <- if(angle_b1 < 0) angle_b1 + 2*pi else angle_b1
  angle_b5 <- if(angle_b5 < 0) angle_b5 + 2*pi else angle_b5
  
  # Calculate the smaller angle between the vectors
  angle_diff <- abs(angle_b5 - angle_b1)
  if(angle_diff > pi) angle_diff <- 2*pi - angle_diff
  
  # Create color vectors for each arrow based on which is minimum
  arrow_colors <- rep("black", 4)
  arrow_colors[min_index] <- "#8FBC8F"
  
  B1 <- min_val
  B1_B5_angle <- angle_diff * 180/pi


  
  # Plot the plane if requested
  if (plot.B1 == T) {
    
    # Create smooth circles for visualization
    circles <- lapply(1:nrow(substi), function(i) {
      # First create basic circle
      circle <- generate_circle(
        center_x = substi$V2[i], 
        center_y = substi$V4[i],
        radius = substi$Radius[i],
        n_points = 100
      )
      
      # Rotate circle points
      rotated <- matrix(nrow = nrow(circle), ncol = 2)
      for (row in 1:nrow(circle)) {
        rotated[row, ] <- aperm(rot.mat %*% as.numeric(circle[row, ]))
      }
      
      data.table::data.table(
        x = rotated[,1],
        y = rotated[,2],
        group = i
      )
    })
    
    # Base plot
    plot <- ggplot2::ggplot(data.frame(tc), ggplot2::aes(x = X1, y = X2)) +
      #ggplot2::geom_point(color = "cadetblue2", shape = 19) +
      ggplot2::geom_vline(xintercept = max_x, color = "darkred", linetype = "dashed") +
      ggplot2::geom_vline(xintercept = min_x, color = "darkred", linetype = "dashed") +
      ggplot2::geom_hline(yintercept = max_y, color = "darkgreen", linetype = "dashed") +
      ggplot2::geom_hline(yintercept = min_y, color = "darkgreen", linetype = "dashed") +
      ggtitle(paste("Plane rotated", i, "degrees"))
    
    # Circles
    plot <- plot +
      ggplot2::geom_polygon(data = do.call(rbind, circles), 
                            ggplot2::aes(x = x, y = y, group = group),
                            fill = "cadetblue2", 
                            color = "cadetblue3",
                            alpha = 0.5) 
    # Add all arrows
    plot <- plot +
      # First arrow (max_x)
      ggplot2::annotate("segment", x = 0, y = 0, xend = max_x, yend = 0, 
                        arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(0.2, "inches")), 
                        color = arrow_colors[1],
                        size = 1.2) +
      # Second arrow (min_x)
      ggplot2::annotate("segment", x = 0, y = 0, xend = min_x, yend = 0, 
                        arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(0.2, "inches")), 
                        color = arrow_colors[2],
                        size = 1.2) +
      # Third arrow (max_y)
      ggplot2::annotate("segment", x = 0, y = 0, xend = 0, yend = max_y, 
                        arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(0.2, "inches")), 
                        color = arrow_colors[3],
                        size = 1.2) +
      # Fourth arrow (min_y)
      ggplot2::annotate("segment", x = 0, y = 0, xend = 0, yend = min_y, 
                        arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(0.2, "inches")), 
                        color = arrow_colors[4],
                        size = 1.2) +
      # B5 arrow
      ggplot2::annotate("segment", x = 0, y = 0, 
                        xend = b5.point[1], 
                        yend = b5.point[2], 
                        arrow = ggplot2::arrow(type = "closed", length = ggplot2::unit(0.2, "inches")), 
                        color = "#CD3333",
                        size = 1.2)
    
    # Add text annotations conditionally
    # Single annotation for the minimum value arrow (B1)
    plot <- plot + switch(min_index,
                          # Case 1: max_x is minimum
                          ggplot2::annotate("text", x = max_x / 2, y = 0, 
                                            label = paste(round(max_x, 3)), 
                                            color = "#8FBC8F", vjust = -1,
                                            size = 4.85, fontface = "bold"),
                          # Case 2: min_x is minimum
                          ggplot2::annotate("text", x = min_x / 2, y = 0, 
                                            label = paste(round(abs(min_x), 2)), 
                                            color = "#8FBC8F", vjust = -1,
                                            size = 4.85, fontface = "bold"),
                          # Case 3: max_y is minimum
                          ggplot2::annotate("text", x = 0, y = max_y / 2, 
                                            label = paste(round(max_y, 2)), 
                                            color = "#8FBC8F", hjust = 1,
                                            size = 4.85, fontface = "bold"),
                          # Case 4: min_y is minimum
                          ggplot2::annotate("text", x = 0, y = min_y / 2, 
                                            label = paste(round(abs(min_y), 2)), 
                                            color = "#8FBC8F", hjust = 1,
                                            size = 4.85, fontface = "bold")
    )
    
    # Add B5 and B1 labels and angle arc
    plot <- plot +
      # B5 label
      ggplot2::geom_label(
        data = data.frame(
          x = (3.3 * b5.point[1] / 5),
          y = (3.3 * b5.point[2] / 5),
          label = paste("B5\n", round(b5.value, 2))
        ),
        aes(x = x, y = y, label = label),
        color = "black",
        fill = "white",
        alpha = 0.8,
        hjust = 0.5,
        size = 4.85,
        fontface = "bold",
        label.padding = unit(0.2, "lines"),
        label.r = unit(0.15, "lines"),
        label.size = 0) +
      # B1 label
      ggplot2::geom_label(
        data = data.frame(
          x = (2.8 * b1_coords[1] / 5),
          y = (2.8 * b1_coords[2] / 5),
          label = paste("B1\n", round(min_val, 2))
        ),
        aes(x = x, y = y, label = label),
        color = "black",
        fill = "white",
        alpha = 0.8,
        hjust = 0.5,
        size = 4.85,
        fontface = "bold",
        label.padding = unit(0.2, "lines"),
        label.r = unit(0.15, "lines"),
        label.size = 0) +
      # Angle arc
      ggplot2::annotate("path",
                        x = 0.5 * cos(seq(min(angle_b1, angle_b5), 
                                          max(angle_b1, angle_b5), 
                                          length.out = 100)),
                        y = 0.5 * sin(seq(min(angle_b1, angle_b5), 
                                          max(angle_b1, angle_b5), 
                                          length.out = 100)),
                        color = "gray40",
                        size = 0.8) +
      # Angle label
      ggplot2::geom_label(
        data = data.frame(
          x = 0.8 * cos((min(angle_b1, angle_b5) + max(angle_b1, angle_b5))/2),
          y = 0.8 * sin((min(angle_b1, angle_b5) + max(angle_b1, angle_b5))/2),
          label = paste(round(angle_diff * (180/pi), 1), "\u00b0")
        ),
        aes(x = x, y = y, label = label),
        color = "black",
        fill = "white",
        alpha = 0.8,
        hjust = 0.5,
        size = 4.85,
        fontface = "bold",
        label.padding = unit(0.2, "lines"),
        label.r = unit(0.15, "lines"),
        label.size = 0) +
      # Theme and scales
      ggplot2::theme_minimal() +
      ggplot2::scale_x_continuous(limits = c(min(tc), max(tc))) +
      ggplot2::scale_y_continuous(limits = c(min(tc), max(tc))) +
      ggplot2::coord_fixed(ratio = 1)
    
    print(plot)
  }
  
  return(data.frame(B1, B1_B5_angle))
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

#' Verloop's sterimol values from xyz files
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
steRimol.xyz <- function(mol, coordinates, CPK = T, only_sub = T, drop = NULL, plot.B1 = T) {
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
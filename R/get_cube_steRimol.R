###### Gaussian .cube files handling ######


#' Extract xyz file from Gausssian cube file
#'
#' This function acts on Gaussian cube file, generating an xyz file
#' of the optimized structure
#' Input is the path to a .cube file
#'
#' @param cube_file path to cube file
#' @return write an xyz file with the same name
#' @export
xyz.from.cube <- function(cube_file) {
  unlink(paste(stringr::str_remove(tools::file_path_sans_ext(cube_file),
                                   '.cube'), ".xyz", sep = ""))
  full.cube <- data.frame(data.table::fread(cube_file, fill = T))
  full.cube <- full.cube[-(1:2),]
  n.atoms <- as.numeric(full.cube[1,1])
  structure <- full.cube[1:(4 + n.atoms), 1:5]
  xyz <- structure[5:(4 + n.atoms), c(1, 3:5)]
  xyz[, 2:4] <- apply(xyz[,2:4], 2, as.numeric)
  xyz[, 2:4] <- xyz[, 2:4] * 0.529177249
  names(xyz) <- c('atom', 'x', 'y', 'z')
  suppressMessages(xyz$atom <- plyr::mapvalues(xyz$atom,
                                               from = atomic_symbols$V1,
                                               to = atomic_symbols$V2
  ))
  num.atoms <- nrow(xyz)
  m <- as.data.frame(matrix(NA, ncol = 4, nrow = 2))
  m[1, 1] <- num.atoms
  m[is.na(m)] <- ""
  names(m) <- names(xyz)
  xyz <- rbind(m, xyz)
  new_xyz <- knitr::kable(xyz, format = "simple", row.names = F,
                          col.names = NULL, align = "l")
  new_xyz <- new_xyz[-1]
  new_xyz <- new_xyz[-length(new_xyz)]
  write(new_xyz, paste(stringr::str_remove(tools::file_path_sans_ext(cube_file),
                                           '.cube'), ".xyz", sep = ""))
}

#' Compute cube based sterimol parameters
#'
#' This function acts on Gaussian cube file
#' NOTE - the function is meant to work in a single molecule
#' data folder and does not allow the choice of a cube file to work with.
#' This is by design and is due to the workflow approach of moleculaR.
#'
#' @param cubefile a .cube file
#' @param coordinates primary axis, a two atom character
#' @param only.sub if TRUE (default) will account only for atoms directly bonded
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @param plot if TRUE (default) will generate a 3D plot of molecule, with the
#' density surface and sterimol arrows.
#' @param degrees define the rotational scan, defaults to 90 (same as 360)
#' @param isovalue what minimal density isovalue should be considered, defaults to 0.036
#' @param plot.B1 should a plot of the projected substituent plane (perpendicular to L) and the B1 feature be presented
#' @return sterimol values for a single molecule
#' @aliases cube.sterimol
#' @export
steRimol.cube <- function(cubefile, coordinates, only.sub = T, drop = NULL, plot = T, degrees = 90,
                          isovalue = 0.036, plot.B1 = T) {
  B.scanner <- function(i, plot.b1 = FALSE, full.topo = FALSE, plot.reduced = FALSE) {
    angle <- i * (pi / 180)
    rot.mat <- matrix(c(cos(angle), -sin(angle), 
                        sin(angle), cos(angle)), 2, 2)
    
    tc_full <- t(rot.mat %*% t(as.matrix(plane))) * 0.529177249
    
    if (!full.topo) {
      x_range <- range(tc_full[,1])
      y_range <- range(tc_full[,2])
      margin <- 0.8
      edge_mask <- tc_full[,1] >= (x_range[2] - margin) | 
        tc_full[,1] <= (x_range[1] + margin) |
        tc_full[,2] >= (y_range[2] - margin) | 
        tc_full[,2] <= (y_range[1] + margin)
      tc <- tc_full[edge_mask, ]
    } else {
      tc <- tc_full
    }
    
    # Calculate metrics using tc
    max_x <- max(tc[, 1])
    min_x <- min(tc[, 1])
    max_y <- max(tc[, 2])
    min_y <- min(tc[, 2])
    
    avs <- abs(c(max_x, min_x, max_y, min_y))
    min_val <- min(avs)
    min_index <- which.min(avs)
    
    b1_coords <- switch(min_index,
                        c(max_x, 0),    
                        c(min_x, 0),    
                        c(0, max_y),    
                        c(0, min_y))    
    
    mag <- function(vector) sqrt(sum((vector)^2))
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
    
    b1 <- min_val
    
    if (plot.b1) {
      plot_data <- if(plot.reduced) tc else tc_full
      
      arrow_colors <- rep("black", 4)
      arrow_colors[min_index] <- "#8FBC8F"
      
      plot <- ggplot2::ggplot(data.frame(plot_data), ggplot2::aes(x = X1, y = X2)) +
        ggplot2::geom_point(color = "cadetblue2", shape = 19) +
        ggplot2::geom_vline(xintercept = max_x, color = "darkred", linetype = "dashed") +
        ggplot2::geom_vline(xintercept = min_x, color = "darkred", linetype = "dashed") +
        ggplot2::geom_hline(yintercept = max_y, color = "darkgreen", linetype = "dashed") +
        ggplot2::geom_hline(yintercept = min_y, color = "darkgreen", linetype = "dashed") +
        ggplot2::ggtitle(paste("Plane rotated", i, "degrees", 
                      if(plot.reduced) " (reduced set)" else ""))
      
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
          ggplot2::aes(x = x, y = y, label = label),
          color = "black",
          fill = "white",
          alpha = 0.8,
          hjust = 0.5,
          size = 4.85,
          fontface = "bold",
          label.padding = ggplot2::unit(0.2, "lines"),
          label.r = ggplot2::unit(0.15, "lines"),
          label.size = 0) +
        # B1 label
        ggplot2::geom_label(
          data = data.frame(
            x = (2.8 * b1_coords[1] / 5),
            y = (2.8 * b1_coords[2] / 5),
            label = paste("B1\n", round(min_val, 2))
          ),
          ggplot2::aes(x = x, y = y, label = label),
          color = "black",
          fill = "white",
          alpha = 0.8,
          hjust = 0.5,
          size = 4.85,
          fontface = "bold",
          label.padding = ggplot2::unit(0.2, "lines"),
          label.r = ggplot2::unit(0.15, "lines"),
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
          ggplot2::aes(x = x, y = y, label = label),
          color = "black",
          fill = "white",
          alpha = 0.8,
          hjust = 0.5,
          size = 4.85,
          fontface = "bold",
          label.padding = ggplot2::unit(0.2, "lines"),
          label.r = ggplot2::unit(0.15, "lines"),
          label.size = 0) +
        # Theme and scales
        ggplot2::theme_minimal() +
        ggplot2::scale_x_continuous(limits = c(min(tc), max(tc))) +
        ggplot2::scale_y_continuous(limits = c(min(tc), max(tc))) +
        ggplot2::coord_fixed(ratio = 1)
      
      print(plot)
    }
    
    return(data.frame(b1 = min_val, b5.value = b5.value, b1_b5_angle = angle_diff, b5.point))
  }
  
  dens.to.pt <- function(point) {
    num.pt <- dense.points[point, ]
    return(c(x.origin + num.pt[[1]] * x.size,
             y.origin + num.pt[[2]] * y.size,
             z.origin + num.pt[[3]] * z.size))
  }
  close.atom <- function(dens.pt) {
    which.min(apply(cbind(trans.co[dens.pt, 1] - trans.co[1:n.atoms, 1],
                          trans.co[dens.pt, 2] - trans.co[1:n.atoms, 2],
                          trans.co[dens.pt, 3] - trans.co[1:n.atoms, 3]),
                    MARGIN = 1,
                    FUN = mag))
  }
  vec.organizer <- function(block) {
    vec <- c(t(density[blocks[block]:(blocks[block] + (ceiling(z.steps/6) - 1)), ]))
    if (length(vec) < ceiling(pracma::nthroot(x.steps * y.steps * z.steps, 3))) {
      vec <- append(vec, rep(0, ceiling(pracma::nthroot(x.steps * y.steps * z.steps, 3)) - length(vec)))
    }
    vec
  }
  pt.space.block <- function(x.block) {
    cbind(x.block, which(x.blocks[[x.block]] != 0, arr.ind = T))
  }
  pt.space.block.binder <- function(list) {
    do.call(rbind, lapply(list, pt.space.block))
  }
  name_of_axis <- stringr::str_replace(coordinates, ' ', '_')
  xyz.from.cube(cubefile)
  origin <- as.numeric(unlist(strsplit(coordinates, " "))[[1]])
  direction <- as.numeric(unlist(strsplit(coordinates, " "))[[2]])
  bonds <- extract.connectivity(list.files(pattern = '.xyz'), threshold_distance = 2.2, keep.HB = F)
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
  
  if (only.sub == T) {
    all_paths <- find_paths_with_nodes(bonds,
                                       as.numeric(unlist(stringr::str_split(coordinates, ' ')))[1],
                                       as.numeric(unlist(stringr::str_split(coordinates, ' ')))[2])
    rlev <- unique(unlist(all_paths))
    if (!is.null(drop)) {
      rlev <- rlev[!rlev %in% drop]
    }
  } else {
    all_paths <- find_paths_with_nodes(bonds,
                                       min(bonds),
                                       max(bonds))
    rlev <- unique(unlist(all_paths))
  }
  
  full.cube <- data.frame(data.table::fread(cubefile,
                                            fill = T))
  full.cube <- full.cube[-(1:2),]
  full.cube <- data.frame(apply(full.cube, 2, as.numeric))
  n.atoms <- as.numeric(full.cube[1,1])
  structure <- full.cube[1:(4 + n.atoms), 1:5]
  density <- full.cube[-c(1:(4 + n.atoms)),]
  x.origin <- structure[1, 2]
  y.origin <- structure[1, 3]
  z.origin <- structure[1, 4]
  x.steps <- structure[2, 1]
  y.steps <- structure[3, 1]
  z.steps <- structure[4, 1]
  x.size <- structure[2, 2]
  y.size <- structure[3, 3]
  z.size <- structure[4, 4]
  cube.volume <- x.size*y.size*z.size*(0.529177249)^3
  xyz <- structure[5:(4 + n.atoms), c(1, 3:5)]
  atomic.numbers <- xyz$V1
  names(xyz) <- c('atom', 'x', 'y', 'z')
  suppressMessages(xyz$atom <- plyr::mapvalues(xyz$atom,
                                               from = atomic_symbols$V1,
                                               to = atomic_symbols$V2
  ))
  
  blocks <- seq(1, nrow(density), ceiling(z.steps/6))
  density[is.na(density)] <- 0
  row.dens <- data.frame(do.call(rbind, lapply(1:length(blocks),
                                               function(block) vec.organizer(block))))
  ## cubes with electrons in them
  elec.cubes <- row.dens[row.dens >= isovalue]
  
  ## molecule volume
  mol.vol <- length(elec.cubes)*cube.volume
  
  row.dens[row.dens > isovalue*1.1 | row.dens < isovalue] <- 0
  
  x.blocks <- list()
  for (i in seq(1, nrow(row.dens), y.steps)) {
    x.blocks[[i]] <- row.dens[i:(i + (y.steps - 1)), ]
  }
  x.blocks <- x.blocks[seq(1, nrow(row.dens), y.steps)]
  non.zero.regions <- which(lapply(x.blocks, function(x) sum(colSums(x))) != 0)
  
  dense.points <- pt.space.block.binder(non.zero.regions)
  
  coordinate.space <- data.frame(do.call(rbind, lapply(1:nrow(dense.points), function(x) dens.to.pt(x))))
  names(coordinate.space) <- names(xyz)[2:4]
  coordinate.space <- rbind(xyz[,2:4], coordinate.space)
  
  
  atoms <- strsplit(coordinates, " ")
  unlisted.atoms <- unlist(atoms)
  numeric.atoms <- as.numeric(unlisted.atoms)
  xyz <- coordinate.space
  count.0 <- function(x) sum(cumsum(x != 0) == 0)
  tag <- NULL
  if (count.0(colSums(xyz[1:n.atoms, 1:3])) >= 2) {
    col.1 <- count.0(xyz[,1])
    col.2 <- count.0(xyz[,2])
    col.3 <- count.0(xyz[,3])
    place.num <- which.max(c(col.1, col.2, col.3))
    new.row <- data.frame(matrix(nrow = 1, ncol = 3))
    new.row[1, place.num] <- 1
    new.row[1, -c(place.num)] <- 0
    names(new.row) <- names(xyz)
    xyz <- data.frame(rbind(xyz[1:n.atoms,], new.row, xyz[n.atoms:nrow(xyz),]))
    numeric.atoms[3] <- n.atoms + 1
    tag <- 1
  }
  mag <- function(vector) {
    sqrt(sum((vector)^2))
  }
  new_origin <- xyz[numeric.atoms[[1]], 1:3]
  new_y <- as.numeric((xyz[numeric.atoms[[2]], 1:3] - new_origin) /
                        mag(xyz[numeric.atoms[[2]], 1:3] - new_origin))
  coplane <- as.numeric((xyz[numeric.atoms[[3]], 1:3] - new_origin) /
                          mag(xyz[numeric.atoms[[3]], 1:3] - new_origin))
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
  new_origin <- xyz[numeric.atoms[[1]], 1:3]
  new_coordinates <- matrix(nrow = dim(xyz)[[1]], ncol = 4)
  trans.co <- data.frame(matrix(nrow = dim(xyz)[[1]], ncol = 4))
  for (i in 1:dim(xyz)[[1]]) {
    new_coordinates[i, 1:3] <- as.numeric(xyz[i, 1:3] - new_origin)
    trans.co[i, 1:3] <- aperm(new_basis %*% new_coordinates[i, 1:3])
  }
  if (!is.null(tag)) trans.co <- trans.co[-(n.atoms + 1),]
  trans.co[(n.atoms + 1):nrow(trans.co), 4] <- as.numeric(unlist(lapply((n.atoms + 1):nrow(trans.co),
                                                                        function(x) close.atom(x))))
  
  trans.co <- data.frame(round(trans.co, 6))
  tcs <- trans.co[trans.co$X4 %in% rlev,1:3]
  tcs.rest.of.mol <- trans.co[-(trans.co$X4 %in% rlev),1:3]
  
  if (plot == T) rgl::plot3d(tcs, col = 'black', axes = F) # Plot for maintenance
  if (plot == T) rgl::plot3d(tcs.rest.of.mol, col = 'grey', add = T)
  if (plot == T) rgl::plot3d(trans.co[1:n.atoms,1:3], col = col_and_size$V2[atomic.numbers],
                             size = col_and_size$V3[atomic.numbers], type = 's', add = T)
  plot.edges <- unlist(lapply(1:nrow(bonds), function(i) as.list(bonds[i,])))
  if (plot == T) rgl::segments3d(trans.co[plot.edges,1:3], add = T, lwd = 3.5)
  if (plot == T) rgl::aspect3d('iso')
  mag <- function(vector) {
    sqrt(vector[[1]]^2 + vector[[2]]^2)
  }
  tcs <- dplyr::mutate(tcs, B = (mag(tcs[,c(1,3)]) * 0.529177249^2))
  L <- max(tcs[,2]) * 0.529177249
  plane <- tcs[, c(1, 3)]
  b1s <- do.call(rbind, lapply(1:degrees, function(row) B.scanner(row)))
  B1 <- b1s$b1[which.min(b1s$b1)]
  b1_b5_angle <- b1s$b1_b5_angle[which.min(b1s$b1)] * 180/pi
  B5 <- b1s$b5.value[which.max(b1s$b5.value)]
  if (plot.B1) invisible(B.scanner(which.min(b1s$b1), plot.b1 = T, T, F))
  loc.B5 <- tcs[which.max(tcs[,4]), 2] * 0.529177249
  if (plot == T) rgl::arrow3d(p0 = c(0, loc.B5/0.529177249, 0),
                              p1 = c(tcs[which.max(tcs$B), 1],
                                     loc.B5/0.529177249,
                                     tcs[which.max(tcs$B), 3]), barblen = 0.03,
                              type = "rotation",  col = "darkgreen") # plot B5 vector
  
  if (plot == T) rgl::arrow3d(p0 = c(0, min(tcs[,2]), 0),
                              p1 = c(0,
                                     max(tcs[,2]),
                                     0), barblen = 0.03,
                              type = "rotation",  col = "darkblue") # plot L vector
  result <- data.frame(round(data.frame(B1, B5, L, loc.B5,b1_b5_angle, mol.vol), 3))
  unlink(paste0(tools::file_path_sans_ext(cubefile), '.xyz'))
  names(result) <- paste0(names(result), '_', name_of_axis)
  return(result)
}

#' Compute cube based sterimol parameters
#'
#' The function expects to be executed in a working directory
#' with cube files. It is recommended
#'
#' @param coor.atoms primary axis, a two atom character
#' @param only.sub if TRUE (default) will account only for atoms directly bonded
#' @param drop numeric value of an atom, from which onward atoms will be dropped from calculation
#' @param plot if TRUE (default) will generate a 3D plot of molecule, with the
#' density surface and sterimol arrows.
#' @param degrees define the rotational scan, defaults to 90 (same as 360)
#' @param isovalue what minimal density isovalue should be considered, defaults to 0.003
#' @return sterimol values for a single molecule
#' @export
steRimol.cube.df <- function(coor.atoms,
                             only.sub = T,
                             drop = NULL,
                             plot = F,
                             degrees = 90,
                             isovalue = 0.003) {
  unlink(list.files(pattern = '.xyz'))
  molecules <- list.files(pattern = '.cube')
  steri.list <- list()
  for (molecule in molecules) {
    steri.list[[match(molecule, molecules)]] <- steRimol.cube(molecule,
                                                              coor.atoms,
                                                              only.sub,
                                                              drop, plot,
                                                              degrees,
                                                              isovalue,
                                                              F)
  }
  steri.dafr <- data.frame(data.table::rbindlist(steri.list))
  row.names(steri.dafr) <- stringr::str_remove(molecules, '.cube')
  return(steri.dafr)
}

#' Calculate cube based sterimol parameters for multiple coordinate sets
#'
#' @param coordinates_vector Vector of coordinate strings (e.g. c("1 2", "3 4"))
#' @param only.sub If TRUE (default) will account only for atoms directly bonded
#' @param drop Numeric value of atoms to drop from calculation
#' @param degrees Define rotational scan, defaults to 90
#' @param isovalue Minimum density isovalue to consider, defaults to 0.036
#' @return Data frame with sterimol values for each coordinate set
#' @export
steRimol.cube.multi <- function(coordinates_vector,
                                only.sub = T,
                                drop = NULL,
                                degrees = 90,
                                isovalue = 0.036) {
  
  multi.df <- lapply(coordinates_vector,
                     function(x) {
                       steRimol.cube.df(x,
                                        only.sub,
                                        drop,
                                        degrees,
                                        isovalue)
                     }
  )
  
  # Combine results and remove duplicate columns
  multi.df.result <- data.frame(do.call(cbind, multi.df))
  multi.df.result <- multi.df.result[!duplicated(as.list(multi.df.result))]
  
  return(multi.df.result)
}


####### ----------------------------------------------------#####
####### -----------------Utility Functions------------------#####
####### ----------------------------------------------------#####

# find_circular_paths - finds 6 membered rings for ring vibrations
#' find_circular_paths - find 6 membered circular paths that include the 
#' primary ring_vib atom. Used in the process of ring vibrations calculation
#' @param bonds bonds dataframe
#' @param node primary ring_vib atom
#' @keywords internal
#' @return 6-membered circular paths that include the primary atom
find_circular_paths <- function(bonds, node) {
  
  # Convert the dataframe to an adjacency matrix
  graph <- matrix(0, nrow = max(bonds), ncol = max(bonds))
  for (i in 1:nrow(bonds)) {
    graph[bonds[i, 1], bonds[i, 2]] <- 1
    graph[bonds[i, 2], bonds[i, 1]] <- 1
  }
  
  first_cycle <- NULL
  
  # Depth-first search function to find cycles
  dfs <- function(start_node, current_node, current_path) {
    current_path <- c(current_path, current_node)
    
    # If the current path forms a cycle of length 7, save it and stop searching
    if (length(current_path) == 7 && current_node == start_node) { # Length 8 because it includes the start node twice
      first_cycle <<- current_path
      return(TRUE)
    }
    
    # Recursive DFS for all neighbors
    neighbors <- which(graph[current_node, ] == 1)
    for (neighbor in neighbors) {
      if (!(neighbor %in% current_path[-1])) {  # Avoid revisiting nodes except the start node
        if (dfs(start_node, neighbor, current_path)) {
          return(TRUE)
        }
      }
    }
    
    return(FALSE)
  }
  
  # Start DFS from the specified node
  dfs(node, node, numeric(0))
  
  if (!is.null(first_cycle)) {
    # Reorder the nodes in the specified order
    reordered_cycle <- c(
      first_cycle[4],   # Fourth node first
      first_cycle[1],   # First node second
      first_cycle[2],   # Second node remains in place
      first_cycle[6],   # Sixth node third
      first_cycle[3],   # Third node fifth
      first_cycle[5]    # Fifth node sixth
    )
    
    # Convert to a character vector with space-separated nodes
    result <- paste(reordered_cycle, collapse = " ")
    return(result)
  } else {
    return(NULL)
  }
}

# cross_product - well... returns the CP result vector
#' compute the cross product of two vetors
#' @param v1 vector 1
#' @param v2 vector 2
#' @keywords internal
#' @return cross product vector
cross_product <- function(v1, v2) {
  as.matrix(as.numeric(c(
    v1[2] * v2[3] - v1[3] * v2[2],
    v1[3] * v2[1] - v1[1] * v2[3],
    v1[1] * v2[2] - v1[2] * v2[1]
  )))
}

# find.center.atom - support for bending vibration extraction, verifies that two
#   atoms share a bond with a center atom, before calculating their bending
#   frequency. Uses an adjency list created from connectivity.

#' Find the atom connecting two atoms, if there is one
#' @param atom1 atom 1 (numeric)
#' @param atom2 atom 2 (numeric)
#' @keywords internal
#' @return center atom or FALSE
find.center.atom <- function(atom1, atom2) {
  bond_pairs <- extract.connectivity()
  adjacency_list <- create_adjacency_list(bond_pairs)
  queue <- c(atom1, atom2)
  visited <- rep(FALSE, length(adjacency_list))
  while (length(queue) > 0) {
    atom <- queue[1]
    queue <- queue[-1]
    if (!visited[atom]) {
      visited[atom] <- TRUE
      for (neighbor in adjacency_list[[atom]]) {
        if (!visited[neighbor]) {
          queue <- c(queue, neighbor)
          if (check.bond.exist(atom1, neighbor) & check.bond.exist(atom2, neighbor)) {
            return(TRUE)
          }
        }
      }
    }
  }
  return(FALSE)
}

#' Create adjacency matrix
#' @param bond_pairs a two column data frame of bonded atoms
#' @keywords internal
#' @return distance matrix for all atoms in a molecule
create_adjacency_list <- function(bond_pairs) {
  adjacency_list <- vector("list", max(bond_pairs))
  for (i in 1:nrow(bond_pairs)) {
    atom1 <- bond_pairs[i,1]
    atom2 <- bond_pairs[i,2]
    adjacency_list[[atom1]] <- c(adjacency_list[[atom1]], atom2)
    adjacency_list[[atom2]] <- c(adjacency_list[[atom2]], atom1)
  }
  return(adjacency_list)
}

# find.bend.freq - support for bending vibration extraction, extarcts a single
#   bending frequency for a given atom pair (assuming they share a center).
#   True functionality comes in the form of several pairs and that of entire sets
#   as user functions. This serves as a utility function, not to be used
#   as a standalone.

#' Find a single bending frequency for a pair of atoms
#' @param atom_pair a character of two atoms (e.g. '1 2')
#' @param threshold minimum frequency value (above fingerprint region)
#' @keywords internal
#' @return frequency
find.bend.freq <- function(atom_pair, threshold = 1350) {
  center.atom <- find.center.atom(atom_pair[1], atom_pair[2])
  if (!center.atom) {
    return('Atoms do not share a center')
  } else {
    mag <- function(vector) {
      sqrt(vector[[1]]^2 + vector[[2]]^2 + vector[[3]]^2)
    }
    vib.info <- data.frame(data.table::fread(list.files(pattern = 'IR'),
                                             col.names = c('Frquency',
                                                           'Intensity'),
                                             colClasses = rep('numeric', 2)))
    atom.1.vecs <- atom.vectors(
      list.files(pattern = paste0('_', as.character(atom_pair[1]), '\\.csv$')))
    atom.2.vecs <- atom.vectors(
      list.files(pattern = paste0('_', as.character(atom_pair[2]), '\\.csv$')))
    cp.mag <- function(row) {
      mag(t(cross_product(atom.1.vecs[row,], atom.2.vecs[row,])))
    }
    freq.num <- which.max(lapply(as.vector(which(vib.info[, 1] > threshold & vib.info[, 1] < 3000)),
                                 cp.mag))
    freq.max <- as.vector(which(vib.info[, 1] > threshold & vib.info[, 1] < 3000))[freq.num]
    return(vib.info[freq.max, 1])
  }
}

# name.vib - helper function for naming
#' helper function for naming
#' @param x bond number
#' @param atom_pairs bonds
#' @keywords internal
#' @return names
name.vib <- function(x, atom_pairs) {
  sub("(\\d+)-(\\d+)","\\1V\\2",
      gsub("\\D+", "-", as.character(atom.pairs.num(atom_pairs)[x])))
}

# check.bond.exist - runs in molecular folder
#   input - comma separated two atom indices. Returns T/F.

#' verifies that a bond exists
#' @param atom1 atom 1
#' @param atom2 atom 2
#' @keywords internal
#' @return logical
check.bond.exist <- function(atom1, atom2){
  bond_pairs <- extract.connectivity()
  for (i in 1:nrow(bond_pairs)) {
    if ((bond_pairs[[i,1]] == atom1 & bond_pairs[[i,2]] == atom2) |
        (bond_pairs[[i,1]] == atom2 & bond_pairs[[i,2]] == atom1)) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# atom.vectors - retrives movement vectors for an atom,
#   based on its index (name of file is always vib_<index>.csv)
#' Retrieve movement vectors for an atom
#' @param file a .feather file
#' @keywords internal
#' @return data frame of movement vectors for a single atom
atom.vectors <- function(file) {
  data.frame(
    data.table::fread(file,
                      col.names = c('x', 'y', 'z'),
                      colClasses = rep('numeric', 3)))
}

####### ----------------------------------------------------#####
####### -------------------User Functions-------------------#####
####### ----------------------------------------------------#####


# stretch.vib and stretch.vib.df - identify bond vibrational frequencies by
#   atom pairs. Expects input of atom.pairs; a character of all
#   atom pairs, by order and separated with a space.
#   example input: '1 3 2 18 27 29' will compute frequencies
#   for 1-2, 2-18 and 27-29.

#' Get a bond's stretching frequency
#'
#' @param atom_pairs a character of atom pairs
#' @param threshold minimal frequency value (above fingerprint region)
#' @keywords internal
#' @return A data frame with stretching frequencies
stretch.vib <- function(atom_pairs, threshold = 1350) {
  xyz <- data.frame(data.table::fread(list.files(pattern = '.xyz'),
                                      col.names = c('atom', 'x', 'y', 'z')))
  xyz <- xyz[, -1]
  xyz <- data.frame(sapply(xyz, as.numeric))
  files <- list.files(pattern = "vib")
  paired <- atom.pairs.num(atom_pairs)
  if (all(unlist(lapply(1:length(paired),
                        function(x) check.bond.exist(paired[[x]][1],
                                                     paired[[x]][2]))))) {
    vec.list <- list()
    for (file in files) {
      vec.list[[match(file, files)]] <- atom.vectors(file)
      names(vec.list)[[match(file, files)]] <- gsub("\\D+", "",
                                                    stringr::str_sub(file,
                                                                     start = 5,
                                                                     end = 7))
    }
    units <- list()
    for (i in 1:length(paired)) {
      units[[i]] <- cbind(
        vec.list[[as.character(paired[[i]][[1]])]],
        vec.list[[as.character(paired[[i]][[2]])]]
      )
      units[[i]] <- data.frame(units[[i]])
      units[[i]][, 7:9] <- xyz[paired[[i]][[1]], ] - xyz[paired[[i]][[2]], ]
    }
    info <- data.frame(data.table::fread(list.files(pattern = 'IR'),
                                         col.names = c('Frquency',
                                                       'Intensity'),
                                         colClasses = rep('numeric', 2)))
    max.list <- list()
    for (i in 1:length(units)) {
      units[[i]][, 10] <- NA
      units[[i]][, 11] <- info[, 1]
      units[[i]] <- units[[i]][units[[i]]$V11 > threshold, ]
      rows <- dim(units[[i]])[1]
      for (j in 1:rows) {
        units[[i]][j, 10] <- abs(pracma::dot(as.matrix(units[[i]][j, 1:3]),
                                             as.matrix(units[[i]][j, 7:9]))) +
          abs(pracma::dot(as.matrix(units[[i]][j, 4:6]),
                          as.matrix(units[[i]][j, 7:9])))
      }
    }
    for (i in 1:length(units)) {
      max.list[[i]] <- which.max(units[[i]][, 10])
    }
    max.list <- data.frame(max.list)
    uni <- data.frame(matrix(nrow = 1, ncol = dim(max.list)[2]))
    for (i in 1:dim(max.list)[2]) {
      uni[1, i] <- units[[i]][max.list[1, i], 11]
    }
    for (i in 1:dim(uni)[2]) {
      colnames(uni)[[i]] <- gsub("\\D+", "-", as.character(paired[i]))
    }
    return(uni)
  } else {
    getwd()
    cat('
One or more bonds do not exist, please validate your choice')
  }
}

#' Get a bond's stretching frequency for a set of molecules
#'
#' @param atom.pairs a character of atom pairs
#' @param threshold minimal wave number value (above fingerprint region)
#' @return A data frame with stretching frequencies
#' @export
stretch.vib.df <- function(atom.pairs, threshold = 1350) {
  molecules <- list.files()
  mol.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    mol.list[[match(molecule, molecules)]] <- stretch.vib(atom.pairs, threshold)
    setwd('..')
  }
  vib.df <- do.call(rbind, mol.list)
  row.names(vib.df) <- molecules
  return(vib.df)
}

# ring.vib, ring.vib.df and ring.vib.df.multi
#   - identify ring characteristic vibrational
#   frequencies by atom pairs.
#   Expects input of atom.pairs; a character of all
#   ring atoms (only 6 membered rings) separated with a space.
#   order matters, assignment does not - as the function relies
#   on atom triplets, it is only imporatnt to maintian a frame of
#   reference. If our primary atom is X, then para will be that on
#   opposite, ortho will be the atoms closest to it and meta, well
#   meta to it (or one atom away from it).
#   ring.vib.df.multi evapluates multiple (or single rings) -
#   this is the offcial function in the package, and sholud be used always
#   when dealing with entire sets. In cases where only one molecule
#   is of interest, use ring.vib.

#' Get a ring's characteristic vibrations frequencies
#'
#' @param primary_atom The atom index of the first ring atom bonded to the basic
#' strcuture.
#' @keywords internal
#' @return A data frame with stretching frequencies
ring.vib <- function(primary_atom) { # single molecule, single ring
  if (is.character(primary_atom)) primary_atom <- as.numeric(primary_atom)
  bonds <- extract.connectivity(list.files(pattern = '.xyz'))
  ordered_ring_atoms <- find_circular_paths(bonds, primary_atom)
  freq <- data.frame(data.table::fread(list.files(pattern = 'IR'),
                                               col.names = c('Frquency',
                                                             'Intensity'),
                                               colClasses = rep('numeric', 2)))
  atoms_info <- lapply(atom.pairs.num(ordered_ring_atoms), as.character)
  xyz <- data.frame(data.table::fread(list.files(pattern = '.xyz'),
                                      col.names = c('atom', 'x', 'y', 'z')))
  xyz <- xyz[, -1]
  xyz <- sapply(xyz, as.numeric)
  xyz <- xyz[unlist(atom.pairs.num(ordered_ring_atoms)), ]
  pa <- xyz[1, 1:3] - xyz[2, 1:3]
  atom.one <- atom.vectors(list.files(pattern = paste0("_", atoms_info[[1]][[1]], '.csv')))
  atom.two <- atom.vectors(list.files(pattern = paste0("_", atoms_info[[3]][[1]], '.csv')))
  atom.three <- atom.vectors(list.files(pattern = paste0("_", atoms_info[[2]][[1]], '.csv')))
  atom.four <- atom.vectors(list.files(pattern = paste0("_", atoms_info[[1]][[2]], '.csv')))
  atom.five <- atom.vectors(list.files(pattern = paste0("_", atoms_info[[2]][[2]], '.csv')))
  atom.six <- atom.vectors(list.files(pattern = paste0("_", atoms_info[[3]][[2]], '.csv')))
  sum.1.3.5 <- list()
  for (i in 1:dim(atom.one)[1]) {
    sum.1.3.5[[i]] <- atom.one[i, ] + atom.three[i, ] + atom.five[i, ]
  }
  vec.sum.1.3.5 <- as.matrix(do.call(rbind, sum.1.3.5))
  sum.2.4.6 <- list()
  for (i in 1:dim(atom.two)[1]) {
    sum.2.4.6[[i]] <- atom.two[i, ] + atom.four[i, ] + atom.six[i, ]
  }
  vec.sum.2.4.6 <- as.matrix(do.call(rbind, sum.2.4.6))
  prods <- list()
  for (i in 1:dim(vec.sum.1.3.5)[1]) {
    prods[[i]] <- pracma::dot(vec.sum.1.3.5[i, ], vec.sum.2.4.6[i, ])
  }
  prod.vec.sums <- data.frame(do.call(rbind, prods))
  prod.vec.sums[, 2] <- freq[, 1]
  for (i in 1:dim(prod.vec.sums)[1]) {
    prod.vec.sums[i, 3] <- abs(sin(angle(vec.sum.2.4.6[i, ], pa)))
  }
  pvec.filter <- dplyr::filter(prod.vec.sums, prod.vec.sums[, 1] != 0)
  vec.prod.filtered <- dplyr::filter(pvec.filter, pvec.filter$V2 > 1550 &
                                       abs(pvec.filter[, 1]) > 0.01)
  dupli.check <- duplicated(vec.prod.filtered$V3) |
    duplicated(vec.prod.filtered$V3, fromLast = T)
  if (any(dupli.check)) {
    right.one <- which.max(abs(vec.prod.filtered[dupli.check,][, 1]))
    unduplicated <- vec.prod.filtered[dupli.check,][right.one,]
    vec.prod.filtered <- rbind(unduplicated, vec.prod.filtered[!dupli.check,])
  }
  vec.prod.filtered[, 1] <- abs(vec.prod.filtered[, 1])
  names(vec.prod.filtered)[1] <- 'dot.prod'
  vec.prod.filtered <- dplyr::arrange(vec.prod.filtered, desc(dot.prod))[1:2, ]
  
  
  result <- data.frame(
    vec.prod.filtered[which.max(vec.prod.filtered[, 3]), 2],
    asin(max(vec.prod.filtered[, 3]))*(180/pi),
    vec.prod.filtered[which.min(vec.prod.filtered[, 3]), 2],
    asin(min(vec.prod.filtered[, 3]))*(180/pi)
  )
  colnames(result) <- c("cross",'cross.angle', "para", 'para.angle')
  if (dim(result)[1] == 0) {
    vec.prod.filtered <- dplyr::filter(pvec.filter, pvec.filter$V2 > 1550 &
                                         pvec.filter$V2 < 1750)
    result <- data.frame(
      vec.prod.filtered[which.max(vec.prod.filtered[, 3]), 2],
      asin(max(vec.prod.filtered[, 3]))*(180/pi),
      vec.prod.filtered[which.min(vec.prod.filtered[, 3]), 2],
      asin(min(vec.prod.filtered[, 3]))*(180/pi)
    )
    colnames(result) <- c("cross",'cross.angle', "para", 'para.angle')
    print(basename(getwd()))
    print('Dot products are lower than 0.01 - returning vibrations in the 1750 - 1550 1/cm region')
  }
  return(result)
}

#' Get a ring's characteristic vibrations frequencies, for a set of molecules
#'
#' @param input The atom index of the first ring atom bonded to the basic
#' strcuture.
#' @return A data frame with stretching frequencies
#' @export
ring.vib.df <- function(input = NA) {
  molecules <- list.dirs(recursive = F, full.names = F)
  if (is.na(input)) input <- readline("Ring atoms - by order -> primary axis (para first), ortho atoms and meta atoms: ")
  ring.vib.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    ring.vib.list[[match(molecule, molecules)]] <- ring.vib(input)
    setwd('..')
  }
  ring.vib.dafr <- data.frame(data.table::rbindlist(ring.vib.list))
  row.names(ring.vib.dafr) <- molecules
  return(ring.vib.dafr)
}

#' Get a ring's characteristic vibrations frequencies
#'
#' used on a set of molecules, applicable for as many rings wanted
#'
#' @param inputs_vector vector of atom indices of the first ring atom bonded to the basic
#' strcuture, for each ring. 
#' 
#' @return A data frame with stretching frequencies
#' @export
ring.vib.multi <- function(inputs_vector) {
  multi.df <- lapply(inputs_vector, ring.vib.df)
  for (i in 1:length(multi.df)) {
    names(multi.df[[i]]) <- paste0(names(multi.df[[i]]), paste0('_', as.character(i)))
  }
  multi.df.result <- data.frame(do.call(cbind, multi.df))
  return(multi.df.result)
}

# bend.vib.set - identify bending frequencies by atom pair.
#   Threshold is set to 1350 on batch uses, it is physically irrelevant to
#   measure vibrational frequencies lower than this threshold.
#   Expects input of atom.pairs; a character of all
#   atom pairs, by order and separated with a space.
#   example input: '1 3 2 18 27 29' will compute frequencies
#   for 1-2, 2-18 and 27-29.

#' Get bending vibrations frequencies
#'
#' @param atom_pairs a character of atom pairs
#' @keywords internal
#' @return A data frame with bending frequencies
bend.vib.set <- function(atom_pairs) {
  result <- data.frame(lapply(atom.pairs.num(atom_pairs), find.bend.freq))
  for (i in 1:ncol(result)) names(result)[i] <-  name.vib(i, atom_pairs)
  return(result)
}

# bond.vib.df - performs bending vibration frequencies for a set of
#   atom pairs, across a molecular library.

#' Get bending vibrations frequencies for a set of molecules
#'
#' @param atom_pairs a character of atom pairs
#' @return A data frame with bending frequencies
#' @export
bend.vib.df <- function(atom_pairs = NA) {
  if (is.na(atom_pairs)) atom_pairs <- readline("Your atom pairs: ")
  molecules <- list.files()
  df <- data.frame(matrix(ncol = length(atom.pairs.num(atom_pairs)),
                          nrow = length(molecules)))
  for (mol in molecules) {
    setwd(mol)
    df[match(mol, molecules), ] <- bend.vib.set(atom_pairs)
    setwd('..')
  }
  for (i in 1:ncol(df)) names(df)[i] <-  name.vib(i, atom_pairs)
  row.names(df) <- molecules
  return(df)
}

#' Get coupled rings characteristic vibrations' frequencies
#'
#' @param primary_atom_1 The atom index of the first ring atom bonded to the basic
#' strcuture.
#' @param primary_atom_2 The atom index of the second ring atom bonded to the basic
#' strcuture.
#' @keywords internal
#' @return A data frame with stretching frequencies
ring.vib.conjugated <- function(primary_atom_1, primary_atom_2) { 
  if (is.character(primary_atom_1)) primary_atom_1 <- as.numeric(primary_atom_1)
  if (is.character(primary_atom_2)) primary_atom_2 <- as.numeric(primary_atom_2)
  bonds <- extract.connectivity(list.files(pattern = '.xyz'))
  ordered_ring_atoms_1 <- find_circular_paths(bonds, primary_atom_1)
  ordered_ring_atoms_2 <- find_circular_paths(bonds, primary_atom_2)
  freq <- data.frame(data.table::fread(list.files(pattern = 'IR'),
                                       col.names = c('Frquency',
                                                     'Intensity'),
                                       colClasses = rep('numeric', 2)))
  atoms_info_1 <- lapply(atom.pairs.num(ordered_ring_atoms_1), as.character)
  atoms_info_2 <- lapply(atom.pairs.num(ordered_ring_atoms_2), as.character)
  
  atom.one.1 <- atom.vectors(list.files(pattern = paste0("_", atoms_info_1[[1]][[1]], '.csv')))
  atom.two.1 <- atom.vectors(list.files(pattern = paste0("_", atoms_info_1[[3]][[1]], '.csv')))
  atom.three.1 <- atom.vectors(list.files(pattern = paste0("_", atoms_info_1[[2]][[1]], '.csv')))
  atom.four.1 <- atom.vectors(list.files(pattern = paste0("_", atoms_info_1[[1]][[2]], '.csv')))
  atom.five.1 <- atom.vectors(list.files(pattern = paste0("_", atoms_info_1[[2]][[2]], '.csv')))
  atom.six.1 <- atom.vectors(list.files(pattern = paste0("_", atoms_info_1[[3]][[2]], '.csv')))
  
  atom.one.2 <- atom.vectors(list.files(pattern = paste0("_", atoms_info_2[[1]][[1]], '.csv')))
  atom.two.2 <- atom.vectors(list.files(pattern = paste0("_", atoms_info_2[[3]][[1]], '.csv')))
  atom.three.2 <- atom.vectors(list.files(pattern = paste0("_", atoms_info_2[[2]][[1]], '.csv')))
  atom.four.2 <- atom.vectors(list.files(pattern = paste0("_", atoms_info_2[[1]][[2]], '.csv')))
  atom.five.2 <- atom.vectors(list.files(pattern = paste0("_", atoms_info_2[[2]][[2]], '.csv')))
  atom.six.2 <- atom.vectors(list.files(pattern = paste0("_", atoms_info_2[[3]][[2]], '.csv')))
  
  mag <- function(vector) {
    sqrt(vector[[1]]^2 + vector[[2]]^2 + vector[[3]]^2)
  }
  
  all.vecs.1 <- cbind(mag(atom.one.1),
                      mag(atom.two.1),
                      mag(atom.three.1),
                      mag(atom.four.1),
                      mag(atom.five.1),
                      mag(atom.six.1))
  
  all.vecs.2 <- cbind(mag(atom.one.2),
                      mag(atom.two.2),
                      mag(atom.three.2),
                      mag(atom.four.2),
                      mag(atom.five.2),
                      mag(atom.six.2))
  
  all.vecs.sum <- rowSums(cbind(all.vecs.1, all.vecs.2))
  vecs.sum.freq <- data.frame(cbind(all.vecs.sum, freq[, 1]))
  vecs.sum.freq <- vecs.sum.freq[vecs.sum.freq$V2 > 1600 & vecs.sum.freq$V2 < 1700, ]
  vec.sum.order <- dplyr::arrange(vecs.sum.freq, V2)
  
  result <- data.frame(
    c.cross.strong = vec.sum.order[1, 2],                    # First row
    c.cross.weak = vec.sum.order[2, 2],                      # Second row
    c.para.weak = vec.sum.order[nrow(vec.sum.order)-1, 2],   # Second to last row
    c.para.strong = vec.sum.order[nrow(vec.sum.order), 2]    # Last row
  )
  
  colnames(result) <- c("cross.c.high",'cross.c.low', "para.c.low", 'para.c.high')
  return(result)
}

#' Get a ring's characteristic vibrations frequencies, for a set of molecules
#'
#' @param primary_atom_1 The atom index of the first ring atom bonded to the basic
#' strcuture.
#' @param primary_atom_2 The atom index of the second ring atom bonded to the basic
#' strcuture.
#' @return A data frame with stretching frequencies
#' @export
ring.vib.conjugated.df <- function(primary_atom_1, primary_atom_2) {
  molecules <- list.dirs(recursive = F, full.names = F)
  ring.vib.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    ring.vib.list[[match(molecule, molecules)]] <- ring.vib.conjugated(primary_atom_1, primary_atom_2)
    setwd('..')
  }
  ring.vib.dafr <- data.frame(data.table::rbindlist(ring.vib.list))
  row.names(ring.vib.dafr) <- molecules
  return(ring.vib.dafr)
}



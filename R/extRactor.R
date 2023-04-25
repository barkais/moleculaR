
######-----------------------------------------------######
######------------------ extRactoR ------------------######
######-----------------------------------------------######
# A Command line tool for the extraction of
#   data from Gaussian log files.
#   The program loops through all log files of a directory,
#   returning a folder with a feather file for each log file.
#   The information in the feather files is all moleculaR needs
#   to extract and calculate molecular features.
#
#   NOTE - This does not include cube files and therefore does not
#   suffice for the cube based stereo-electronic features.
#   To extract these features, please create cube files for all the molecules
#   you wish to analyze.
#   You will need the .chk files turned into .fchk (formchk), which are
#   processed by Gaussian's cubegen program. This of-course requires you
#   to transfer them to your local machine, in the case you choose
#   to work on it.

#' Pull NBO charges from Gaussian log files
#'
#' This function acts on a read Gaussian log file.
#' No input is needed.
#'
#' No Parameters
#' @keywords internal
#' @return A single column data frame with nbo charges
extRact.NBO <- function() {
  pattern1 <- 'Summary of'
  pattern2 <- '====='
  text <- main
  first_index <- max(grep(pattern1, text))
  last_index <- first_index + min(grep(pattern2,
             text[(first_index + 1):length(text)]))
  text <- text[(first_index + 6):(last_index - 1)]
  extract_nbo_values <- function(string) {
    tokens <- strsplit(string, "\\s+")[[1]]
    numerics <- as.numeric(grep("^-?[0-9]+([.][0-9]+)?$", tokens, value = TRUE))
    return(numerics[2])
  }
  result <- data.frame(do.call(rbind, lapply(text, extract_nbo_values)))
  names(result) <- 'NPA'
  return(result)
}

#' Pull the dipole moment vector from Gaussian log files
#'
#' This function acts on a read Gaussian log file.
#' No input is needed.
#'
#' No Parameters
#' @keywords internal
#' @return A 4 column data frame with dipole moment components and magnitude
extRact.Dipole <- function() {
  text <- main
  index <- max(grep("Dipole moment", text))
  tokens <- strsplit(text[index + 1], "\\s+")[[1]]
  numerics <- as.numeric(grep("^-?[0-9]+([.][0-9]+)?$", tokens, value = TRUE))
  result <- data.frame(t(numerics))
  names(result) <- c('x', 'y', 'z', 'Total')
  return(result)
}

#' Pull xyz of the optimized structure from Gaussian log files
#'
#' This function acts on a read Gaussian log file.
#' No input is needed.
#'
#' No Parameters
#' @return A 4 column data frame with atom symbols and coordinates
#' @keywords internal
extRact.xyz <- function() {
  text <- main
  if (length(grep('Standard orientation', text)) == 0) {
    pattern1 <- 'Input orientation'
  } else {
    pattern1 <- 'Standard orientation'
  }
  pattern2 <- '---'
  first_index <- max(grep(pattern1, text))
  last_index <- first_index +
    grep(pattern2,text[(first_index + 1):length(text)])[3]
  text <- text[(first_index + 5):(last_index - 1)]
  extract_xyz_values <- function(string) {
    tokens <- strsplit(string, "\\s+")[[1]]
    numerics <- as.numeric(grep("^-?[0-9]+([.][0-9]+)?$", tokens, value = TRUE))
    return(numerics[c(2, 4:6)])
  }
  result <- data.frame(do.call(rbind, lapply(text, extract_xyz_values)))
  names(result) <- c('Atom', 'x', 'y', 'z')
  result$Atom <- round(result$Atom, 1)
  suppressMessages(result$Atom <- plyr::mapvalues(result$Atom,
    from = atomic_symbols$V1,
    to = atomic_symbols$V2
  ))
  num.atoms <- nrow(result)
  m <- as.data.frame(matrix(NA, ncol = 4, nrow = 2))
  m[1, 1] <- num.atoms
  m[is.na(m)] <- ""
  names(m) <- names(result)
  result <- rbind(m, result)
  return(result)
}

#' Pull IR spectrum from Gaussian log files
#'
#' This function acts on a read Gaussian log file.
#' No input is needed.
#'
#' No Parameters
#' @keywords internal
#' @return A 2 column data frame with IR vibrational frequencies and intensities
extRact.spectrum <- function() {
  text <- main
  pattern1 <- 'Harmonic'
  pattern2 <- 'Thermochemistry'
  first_index <- grep(pattern1, text)
  last_index <- grep(pattern2,text)
  text <- text[first_index:last_index]
  freqs.text <- text[grep('Frequencies', text)]
  intens.text <- text[grep('IR', text)[-1]]
  extract_freq_IR_values <- function(string) {
    tokens <- strsplit(string, "\\s+")[[1]]
    numerics <- as.numeric(grep("^-?[0-9]+([.][0-9]+)?$", tokens, value = TRUE))
    return(numerics)
  }
  Frequencies <- unlist(lapply(freqs.text, extract_freq_IR_values))
  IR_Intensity <- unlist(lapply(intens.text, extract_freq_IR_values))
  spectrum <- data.frame(Frequencies, IR_Intensity)
  return(spectrum)
}

#' Pull isotropic and anisotropic polarizability from Gaussian log files
#'
#' This function acts on a read Gaussian log file.
#' No input is needed.
#'
#' No Parameters
#' @keywords internal
#' @return A 2 column data frame with polarizabilities
extRact.polarizability <- function() {
  text <- main
  index <- max(grep('Dipole polarizability', text))
  text <- text[index:(index + 10)]
  line.aniso <- grep('aniso', text, value = T)
  tokens.an <- strsplit(line.aniso, "\\s+")[[1]]
  tokens.an <- stringr::str_replace_all(tokens.an, 'D', 'E')
  numerics.an <- as.numeric(tokens.an[3],
                         scientific = T)
  line.iso <- grep('iso', text, value = T)
  tokens.iso <- strsplit(line.iso, "\\s+")[[1]]
  tokens.iso <- stringr::str_replace_all(tokens.iso, 'D', 'E')
  numerics.iso <- as.numeric(tokens.iso[3],
                            scientific = T)
  result <- data.frame(t(c(numerics.iso, numerics.an)))
  return(result)
}

#' Pull atom movement vectors from Gaussian log files
#'
#' This function acts on a read Gaussian log file.
#' No input is needed.
#'
#' No Parameters
#' @keywords internal
#' @return A 3 column * number of atoms column data frame with movements vectors
extRact.vectors <- function() {
  text <- main
  pattern1 <- 'Harmonic'
  pattern2 <- 'Thermochemistry'
  first_index <- grep(pattern1, text)
  last_index <- grep(pattern2,text)
  text <- text[first_index:last_index]
  text <- stringr::str_trim(text)
  start.count <- min(grep("Atom ", text, value = F))
  end.count <- start.count +
    (min(grep(" {10}", text[start.count:length(text)])) - 1)
  n.atoms <- length(text[(start.count + 1):(end.count - 1)])
  find_atom_vectors <- function(atom) {
    patt <- paste0('^', as.character(atom), ' ')
    temp <-  text[grep(patt, text, value = F)]
    remove <- grep(" {10}", temp)
    if (length(remove) > 0) temp <- temp[-remove]
    assign(paste0('vib_', as.character(atom)), temp)
  }
  result <- lapply(1:n.atoms, find_atom_vectors)
  extract_detail_values <- function(string) {
    tokens <- strsplit(string, "\\s+")[[1]]
    numerics <- as.numeric(grep("^-?[0-9]+([.][0-9]+)?$", tokens, value = TRUE))
    return(numerics)
  }
  turn_numeric <- function(atom) {
    df <- data.frame(do.call(rbind, lapply(result[[atom]],
                                           extract_detail_values)))
    df <- df[, 3:11]
    names(df) <- c("x", "y", "z", "x", "y", "z", "x", "y", "z")
    iterations <- dim(df)[1]
    datalist <- list()
    for (i in 1:iterations) {
      a <- df[i, 1:3]
      b <- df[i, 4:6]
      c <- df[i, 7:9]
      d <- rbind(a, b, c)
      datalist[[i]] <- assign(
        paste("out", i, sep = "."),
        data.frame(d, row.names = c(1, 2, 3))
      )
    }
    output <- do.call(rbind, datalist)
    row.names(output) <- seq(1, iterations * 3, 1)
    return(output)
  }
  return(data.frame(lapply(1:n.atoms, turn_numeric)))
}

#' Extract and compress all needed information from Gaussian log files.
#'
#' This function acts on a set of Gaussian log files.
#' No input is needed.
#'
#' No Parameters
#' @return A .feather file for each log file
#' @aliases extRactoR
#' @export
extRactoR <- function() {
  main <- character()
  instruction.1 <- svDialogs::dlg_message('
         Please choose the .log files location.
        ', type = 'ok')

  log.files.dir <- svDialogs::dlg_dir()$res

  output.folder.name <- svDialogs::dlg_input('
         How would you like to call the output folder?'
         )$res

  output.folder.path.default <- svDialogs::dlg_message('
         Would you like to save the output folder in the same
         directory as that of the .log files? (Recommended).
        ', type = 'yesno')$res

  if (output.folder.path.default == 'no') {
    instruction.2 <- svDialogs::dlg_message('
         Please choose a location.
        ', type = 'ok')

    output.folder.path <- svDialogs::dlg_dir()$res

    if (output.folder.name == '-') {
      output.folder.name <- paste0(output.folder.path,
                                   '/',
                                   'moleculaR_Files')
      dir.create(output.folder.name)
    } else {
      output.folder.name <- paste0(output.folder.path,
                                   '/',
                                   output.folder.name)
      dir.create(output.folder.name)
    }
  } else {
    if (output.folder.name == '-') {
      output.folder.name <- 'moleculaR_Files'
      dir.create(output.folder.name)
    } else {
      output.folder.name <- paste0(log.files.dir,
                                   '/',
                                   output.folder.name)
      dir.create(output.folder.name)
    }

  }

  setwd(log.files.dir)

  for (file in list.files(pattern = '.log')) {
    main.big <- data.table::fread(file, sep = "?", header = FALSE, quote="")[[1L]]
    scf.dones <- grep('SCF Done', main.big)
    sec.to.last <- scf.dones[length(scf.dones) - 1]
    main <<- main.big[sec.to.last:length(main.big)]
    raw_data <- (list(
      extRact.xyz(),
      extRact.Dipole(),
      extRact.polarizability(),
      extRact.NBO(),
      extRact.spectrum(),
      extRact.vectors()
    ))
    df.result <- data.frame(
      matrix(
        ncol = sum(unlist(lapply(1:6, function(x) ncol(raw_data[[x]])))),
        nrow = max(unlist(lapply(1:6, function(x) nrow(raw_data[[x]]))))))
    df.result[1:nrow(raw_data[[1]]), 1:4] <- raw_data[[1]]
    df.result[1, 5:8] <- raw_data[[2]]
    df.result[1:nrow(raw_data[[3]]), 9:10] <- raw_data[[3]]
    df.result[1:nrow(raw_data[[4]]), 11] <- raw_data[[4]]
    df.result[1:nrow(raw_data[[5]]), 12:13] <- raw_data[[5]]
    df.result[1:nrow(raw_data[[6]]), 14:ncol(df.result)] <-  raw_data[[6]]
    feather::write_feather(df.result,
                   paste0(output.folder.name,
                          '/',
                          tools::file_path_sans_ext(file),
                          '.feather'))
  }
  setwd(output.folder.name)
  cat('
Done!
      ')
  cat(paste0('
Current working directory is set to ',getwd()))
}

#' Extract and compress all needed information from Gaussian log files - no GUI
#'
#' This function acts on a set of Gaussian log files.
#' No input is needed.
#'
#' No Parameters
#' @return A .feather file for each log file
#' @aliases extRactoR
#' @export
extRactoR.auto <- function() {
  main <- character()
  dir.create('Extracted_info')
  for (file in list.files(pattern = '.log')) {
    main.big <- data.table::fread(file,
                                  sep = "?",
                                  header = FALSE,
                                  quote = "")[[1L]]
    scf.dones <- grep('SCF Done', main.big)
    sec.to.last <- scf.dones[length(scf.dones) - 1]
    main <<- main.big[sec.to.last:length(main.big)]
    raw_data <- (list(
      extRact.xyz(),
      extRact.Dipole(),
      extRact.polarizability(),
      extRact.NBO(),
      extRact.spectrum(),
      extRact.vectors()
    ))
    df.result <- data.frame(
      matrix(
        ncol = sum(unlist(lapply(1:6, function(x) ncol(raw_data[[x]])))),
        nrow = max(unlist(lapply(1:6, function(x) nrow(raw_data[[x]]))))))
    df.result[1:nrow(raw_data[[1]]), 1:4] <- raw_data[[1]]
    df.result[1, 5:8] <- raw_data[[2]]
    df.result[1:nrow(raw_data[[3]]), 9:10] <- raw_data[[3]]
    df.result[1:nrow(raw_data[[4]]), 11] <- raw_data[[4]]
    df.result[1:nrow(raw_data[[5]]), 12:13] <- raw_data[[5]]
    df.result[1:nrow(raw_data[[6]]), 14:ncol(df.result)] <-  raw_data[[6]]
    feather::write_feather(df.result,
                           paste0('extracted_info',
                                  '/',
                                  tools::file_path_sans_ext(file),
                                  '.feather'))
  }
  cat('
Done!
      ')
}



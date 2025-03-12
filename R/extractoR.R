
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
  tryCatch({
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
  }, error = function(e) {
    result <- data.frame(NA)
    names(result) <- 'NPA'
    return(result)
  })
}


#' Pull Hirshfeld charges from Gaussian log files
#'
#' This function acts on a read Gaussian log file.
#' No input is needed.
#'
#' No Parameters
#' @keywords internal
#' @return A single column data frame with Hirshfeld charges
extRact.Hirsh <- function() {
  tryCatch({
    pattern1 <- 'Hirshfeld'
    pattern2 <- 'Tot'
    text <- main
    first_index <- grep(pattern1, text)[length(grep(pattern1, text)) - 1]
    last_index <- first_index + min(grep(pattern2,
                                         text[(first_index + 1):length(text)]))
    text <- text[(first_index + 2):(last_index - 1)]
    extract_hirsh_values <- function(string) {
      tokens <- strsplit(string, "\\s+")[[1]]
      numerics <- as.numeric(grep("^-?[0-9]+([.][0-9]+)?$", tokens, value = TRUE))
      return(numerics[2])
    }
    result <- data.frame(do.call(rbind, lapply(text, extract_hirsh_values)))
    names(result) <- 'Hirshfeld'
    return(result)
  }, error = function(e) {
    result <- data.frame(NA)
    names(result) <- 'Hirshfeld'
    return(result)
  })
}


#' Pull CM5 charges from Gaussian log files
#'
#' This function acts on a read Gaussian log file.
#' No input is needed.
#'
#' No Parameters
#' @keywords internal
#' @return A single column data frame with CM5 charges
extRact.CM5 <- function() {
  tryCatch({
    pattern1 <- 'Hirshfeld'
    pattern2 <- 'Tot'
    text <- main
    first_index <- grep(pattern1, text)[length(grep(pattern1, text)) - 1]
    last_index <- first_index + min(grep(pattern2,
                                         text[(first_index + 1):length(text)]))
    text <- text[(first_index + 2):(last_index - 1)]
    extract_cm5_values <- function(string) {
      tokens <- strsplit(string, "\\s+")[[1]]
      numerics <- as.numeric(grep("^-?[0-9]+([.][0-9]+)?$", tokens, value = TRUE))
      return(numerics[7])
    }
    result <- data.frame(do.call(rbind, lapply(text, extract_cm5_values)))
    names(result) <- 'CM5'
    return(result)
  }, error = function(e) {
    result <- data.frame(NA)
    names(result) <- 'CM5'
    return(result)
  })
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
  first_index <- max(grep(pattern1, text))
  last_index <- max(grep(pattern2,text))
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
  suppressWarnings(tryCatch({
    text <- main
    index <- max(grep('Dipole polarizability', text))
    text <- text[index:(index + 10)]
    
    line.aniso <- grep('aniso', text, value = TRUE)
    tokens.an <- strsplit(line.aniso, "\\s+")[[1]]
    tokens.an <- stringr::str_replace_all(tokens.an, 'D', 'E')
    numerics.an <- as.numeric(tokens.an[complete.cases(suppressWarnings(as.numeric(tokens.an)))],
                              scientific = TRUE)[1]
    
    line.iso <- grep('iso', text, value = TRUE)
    tokens.iso <- strsplit(line.iso, "\\s+")[[1]]
    tokens.iso <- stringr::str_replace_all(tokens.iso, 'D', 'E')
    numerics.iso <- as.numeric(tokens.iso[complete.cases(suppressWarnings(as.numeric(tokens.an)))],
                               scientific = TRUE)[1]
    
    result <- data.frame(t(c(numerics.iso, numerics.an)))
    names(result) <- c('Isotropic_Polar', 'Anisotropic_Polar')
    return(result)
    
  }, error = function(e) {
    data.frame(matrix(nrow = 1, ncol = 2))  # Return NA if any error occurs
  }))
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
  first_index <- max(grep(pattern1, text))
  last_index <- max(grep(pattern2,text))
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

#' Extract and compress all needed information from Gaussian log files
#'
#' This function acts on a set of Gaussian log files.
#' No input is needed.
#'
#' No Parameters
#' @return A .feather file for each log file
#' @aliases extRactoR.auto
#' @export
extractoR <- function() {
  main <- character()
  dir.create('Feather_Files')
  for (file in list.files(pattern = '.log')) {
    message("Processing ", file)
    tryCatch(
      {
        main.big <- data.table::fread(file, sep = "?", header = FALSE, quote="")[[1L]]
        scf.dones <- grep('SCF Done', main.big)
        sec.to.last <- scf.dones[length(scf.dones) - 1]
        main <<- main.big[sec.to.last:length(main.big)]
        raw_data <- (list(
          extRact.xyz(),
          extRact.Dipole(),
          extRact.polarizability(),
          extRact.NBO(),
          extRact.Hirsh(),
          extRact.CM5(),
          extRact.spectrum(),
          extRact.vectors()
        ))
        df.result <- data.frame(
          matrix(
            ncol = sum(unlist(lapply(1:8, function(x) ncol(raw_data[[x]])))),
            nrow = max(unlist(lapply(1:8, function(x) nrow(raw_data[[x]]))))))
        df.result[1:nrow(raw_data[[1]]), 1:4] <- raw_data[[1]]
        df.result[1, 5:8] <- raw_data[[2]]
        df.result[1:nrow(raw_data[[3]]), 9:10] <- raw_data[[3]]
        df.result[1:nrow(raw_data[[4]]), 11] <- raw_data[[4]]
        df.result[1:nrow(raw_data[[5]]), 12] <- raw_data[[5]]
        df.result[1:nrow(raw_data[[6]]), 13] <- raw_data[[6]]
        df.result[1:nrow(raw_data[[7]]), 14:15] <- raw_data[[7]]
        df.result[1:nrow(raw_data[[8]]), 16:ncol(df.result)] <-  raw_data[[8]]
        feather::write_feather(df.result,
                               paste0('Feather_Files',
                                      '/',
                                      tools::file_path_sans_ext(file),
                                      '.feather'))
      }, error = function(e){cat(paste0("Error while extracting from ", file, 
                                        ' - check for errors in the log file.'))}
    )
  }
  message('
Done!
      ')
}

#' Extract and compress all needed information from Gaussian log files (parallel version)
#'
#' This function acts on a set of Gaussian log files using parallel processing to improve performance.
#' The function processes multiple files simultaneously using the parallel package.
#'
#' @param num_cores Number of cores to use for parallel processing. Default is 
#'        (detected cores - 1) or 1, whichever is greater.
#' @param output_dir Directory to save the feather files. Default is 'Feather_Files'.
#'
#' @return A .feather file for each log file in the specified output directory
#' @aliases extRactoR.auto
#' @export
extractoR.parallel <- function(num_cores = max(1, parallel::detectCores() - 1), 
                               output_dir = 'Feather_Files') {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Get list of log files
  log_files <- list.files(pattern = '\\.log$')
  
  if (length(log_files) == 0) {
    message("No .log files found in the current directory.")
    return(invisible(NULL))
  }
  
  # Function to process a single log file
  process_log_file <- function(file) {
    tryCatch({
      # Create a local variable for this process to avoid shared state issues
      main <- character()
      
      # Read the log file
      main.big <- data.table::fread(file, sep = "?", header = FALSE, quote="")[[1L]]
      scf.dones <- grep('SCF Done', main.big)
      
      # If there are not enough SCF Done lines, return with error
      if (length(scf.dones) < 2) {
        return(paste0("Error while extracting from ", file, 
                      " - not enough SCF Done lines found"))
      }
      
      # Use the section from the second-to-last SCF Done to the end
      sec.to.last <- scf.dones[length(scf.dones) - 1]
      main <- main.big[sec.to.last:length(main.big)]
      
      # Extract all data components
      raw_data <- list(
        extRact.xyz(main),
        extRact.Dipole(main),
        extRact.polarizability(main),
        extRact.NBO(main),
        extRact.Hirsh(main),
        extRact.CM5(main),
        extRact.spectrum(main),
        extRact.vectors(main)
      )
      
      # Construct the result data frame
      df.result <- data.frame(
        matrix(
          ncol = sum(unlist(lapply(1:8, function(x) ncol(raw_data[[x]])))),
          nrow = max(unlist(lapply(1:8, function(x) nrow(raw_data[[x]]))))
        )
      )
      
      # Fill the result data frame
      df.result[1:nrow(raw_data[[1]]), 1:4] <- raw_data[[1]]
      df.result[1, 5:8] <- raw_data[[2]]
      df.result[1:nrow(raw_data[[3]]), 9:10] <- raw_data[[3]]
      df.result[1:nrow(raw_data[[4]]), 11] <- raw_data[[4]]
      df.result[1:nrow(raw_data[[5]]), 12] <- raw_data[[5]]
      df.result[1:nrow(raw_data[[6]]), 13] <- raw_data[[6]]
      df.result[1:nrow(raw_data[[7]]), 14:15] <- raw_data[[7]]
      df.result[1:nrow(raw_data[[8]]), 16:ncol(df.result)] <- raw_data[[8]]
      
      # Write the result to a feather file
      output_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(file), '.feather'))
      feather::write_feather(df.result, output_file)
      
      return(paste0("Successfully processed ", file))
      
    }, error = function(e) {
      return(paste0("Error while extracting from ", file, 
                    " - ", e$message))
    })
  }
  
  # Modified extractor functions to accept main as an argument rather than using global variable
  
  extRact.NBO <- function(main) {
    tryCatch({
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
    }, error = function(e) {
      result <- data.frame(NA)
      names(result) <- 'NPA'
      return(result)
    })
  }
  
  extRact.Hirsh <- function(main) {
    tryCatch({
      pattern1 <- 'Hirshfeld'
      pattern2 <- 'Tot'
      text <- main
      first_index <- grep(pattern1, text)[length(grep(pattern1, text)) - 1]
      last_index <- first_index + min(grep(pattern2,
                                           text[(first_index + 1):length(text)]))
      text <- text[(first_index + 2):(last_index - 1)]
      extract_hirsh_values <- function(string) {
        tokens <- strsplit(string, "\\s+")[[1]]
        numerics <- as.numeric(grep("^-?[0-9]+([.][0-9]+)?$", tokens, value = TRUE))
        return(numerics[2])
      }
      result <- data.frame(do.call(rbind, lapply(text, extract_hirsh_values)))
      names(result) <- 'Hirshfeld'
      return(result)
    }, error = function(e) {
      result <- data.frame(NA)
      names(result) <- 'Hirshfeld'
      return(result)
    })
  }
  
  extRact.CM5 <- function(main) {
    tryCatch({
      pattern1 <- 'Hirshfeld'
      pattern2 <- 'Tot'
      text <- main
      first_index <- grep(pattern1, text)[length(grep(pattern1, text)) - 1]
      last_index <- first_index + min(grep(pattern2,
                                           text[(first_index + 1):length(text)]))
      text <- text[(first_index + 2):(last_index - 1)]
      extract_cm5_values <- function(string) {
        tokens <- strsplit(string, "\\s+")[[1]]
        numerics <- as.numeric(grep("^-?[0-9]+([.][0-9]+)?$", tokens, value = TRUE))
        return(numerics[7])
      }
      result <- data.frame(do.call(rbind, lapply(text, extract_cm5_values)))
      names(result) <- 'CM5'
      return(result)
    }, error = function(e) {
      result <- data.frame(NA)
      names(result) <- 'CM5'
      return(result)
    })
  }
  
  extRact.Dipole <- function(main) {
    text <- main
    index <- max(grep("Dipole moment", text))
    tokens <- strsplit(text[index + 1], "\\s+")[[1]]
    numerics <- as.numeric(grep("^-?[0-9]+([.][0-9]+)?$", tokens, value = TRUE))
    result <- data.frame(t(numerics))
    names(result) <- c('x', 'y', 'z', 'Total')
    return(result)
  }
  
  extRact.xyz <- function(main) {
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
    
    # Get atomic_symbols from moleculaR namespace 
    atomic_symbols_data <- get("atomic_symbols", envir = asNamespace("moleculaR"))
    
    suppressMessages(result$Atom <- plyr::mapvalues(result$Atom,
                                                    from = atomic_symbols_data$V1,
                                                    to = atomic_symbols_data$V2
    ))
    num.atoms <- nrow(result)
    m <- as.data.frame(matrix(NA, ncol = 4, nrow = 2))
    m[1, 1] <- num.atoms
    m[is.na(m)] <- ""
    names(m) <- names(result)
    result <- rbind(m, result)
    return(result)
  }
  
  extRact.spectrum <- function(main) {
    text <- main
    pattern1 <- 'Harmonic'
    pattern2 <- 'Thermochemistry'
    first_index <- max(grep(pattern1, text))
    last_index <- max(grep(pattern2,text))
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
  
  extRact.polarizability <- function(main) {
    suppressWarnings(tryCatch({
      text <- main
      index <- max(grep('Dipole polarizability', text))
      text <- text[index:(index + 10)]
      
      line.aniso <- grep('aniso', text, value = TRUE)
      tokens.an <- strsplit(line.aniso, "\\s+")[[1]]
      tokens.an <- stringr::str_replace_all(tokens.an, 'D', 'E')
      numerics.an <- as.numeric(tokens.an[complete.cases(suppressWarnings(as.numeric(tokens.an)))],
                                scientific = TRUE)[1]
      
      line.iso <- grep('iso', text, value = TRUE)
      tokens.iso <- strsplit(line.iso, "\\s+")[[1]]
      tokens.iso <- stringr::str_replace_all(tokens.iso, 'D', 'E')
      numerics.iso <- as.numeric(tokens.iso[complete.cases(suppressWarnings(as.numeric(tokens.an)))],
                                 scientific = TRUE)[1]
      
      result <- data.frame(t(c(numerics.iso, numerics.an)))
      names(result) <- c('Isotropic_Polar', 'Anisotropic_Polar')
      return(result)
      
    }, error = function(e) {
      data.frame(matrix(nrow = 1, ncol = 2))  # Return NA if any error occurs
    }))
  }
  
  extRact.vectors <- function(main) {
    text <- main
    pattern1 <- 'Harmonic'
    pattern2 <- 'Thermochemistry'
    first_index <- max(grep(pattern1, text))
    last_index <- max(grep(pattern2,text))
    text <- text[first_index:last_index]
    text <- stringr::str_trim(text)
    start.count <- min(grep("Atom ", text, value = F))
    end.count <- start.count +
      (min(grep(" {10}", text[start.count:length(text)])) - 1)
    n.atoms <- length(text[(start.count + 1):(end.count - 1)])
    
    result_list <- list()
    for (atom in 1:n.atoms) {
      patt <- paste0('^', as.character(atom), ' ')
      temp <- text[grep(patt, text, value = F)]
      remove <- grep(" {10}", temp)
      if (length(remove) > 0) temp <- temp[-remove]
      result_list[[atom]] <- temp
    }
    
    extract_detail_values <- function(string) {
      tokens <- strsplit(string, "\\s+")[[1]]
      numerics <- as.numeric(grep("^-?[0-9]+([.][0-9]+)?$", tokens, value = TRUE))
      return(numerics)
    }
    
    turn_numeric <- function(atom) {
      df <- data.frame(do.call(rbind, lapply(result_list[[atom]],
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
        datalist[[i]] <- data.frame(d, row.names = c(1, 2, 3))
      }
      output <- do.call(rbind, datalist)
      row.names(output) <- seq(1, iterations * 3, 1)
      return(output)
    }
    
    return(data.frame(lapply(1:n.atoms, turn_numeric)))
  }
  
  # Process files in parallel
  message(sprintf("Processing %d files using %d cores", length(log_files), num_cores))
  
  # Use parallel::mclapply for parallel processing
  results <- parallel::mclapply(
    log_files,
    process_log_file,
    mc.cores = num_cores
  )
  
  # Count successes and errors
  success_count <- sum(grepl("^Successfully", unlist(results)))
  error_count <- sum(grepl("^Error", unlist(results)))
  
  # Display results
  message(sprintf("\nProcessing complete: %d files processed successfully, %d failed", 
                  success_count, error_count))
  
  # Display error messages if any
  if (error_count > 0) {
    message("\nErrors:")
    error_msgs <- results[grepl("^Error", unlist(results))]
    for (msg in error_msgs) {
      message(msg)
    }
  }
  
  message("\nDone!")
  
  invisible(results)
}


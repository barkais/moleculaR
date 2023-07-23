######-----------------------------------------------######
######------------------ unwRapper ------------------######
######-----------------------------------------------######
#
# A base function designed to unwrap feather files created by
#   extRactoR.

#' Single molecule unwrapper for feather files
#'
#' Translate stored information into human readable csv files
#' @param feather_file A feather file resulting from extRactoR
#' @keywords internal
#' @return folder with optimized structures as xyz, and information needed for moleculaR.
unwRapper.single <- function(feather_file) {
  data <- feather::read_feather(feather_file)
  xyz <- data[complete.cases(data[, 1:4]), 1:4]
  dipole <- data[complete.cases(data[, 5:8]), 5:8]
  polar <- data.frame(data[complete.cases(data[, 9:10]), 9:10])
  nbo <- data.frame(data[complete.cases(data[, 11]), 11])
  ir.spectrum <- data[complete.cases(data[, 12:13]), 12:13]
  vectors <- data[complete.cases(data[, 14:ncol(data)]), 14:ncol(data)]

  # Define naming for all files
  name <- tools::file_path_sans_ext(feather_file)

  if (dir.exists(name)) {
    # If the folder exists, delete it
    unlink(name, recursive = TRUE)
  }

  dir.create(name) # create a folder for the csv files

  # write xyz file
  new_xyz <- knitr::kable(xyz,
                          format = "simple",
                          row.names = F,
                          col.names = NULL,
                          align = "l")
  new_xyz <- new_xyz[-1]
  new_xyz <- new_xyz[-length(new_xyz)]
  write(new_xyz, paste0(name, "/opt_structure.xyz"))

  # write dipole csv
  write.csv(dipole, paste0(name, '/dipole.csv'),
            row.names = FALSE)

  # write polarizability csv
  write.csv(polar, paste0(name, '/polar.csv'),
            row.names = FALSE)

  # write NBO csv
  write.csv(nbo, paste0(name, '/nbo.csv'),
            row.names = FALSE)

  # write IR spectrum csv
  write.csv(ir.spectrum, paste0(name, '/IR.csv'),
            row.names = FALSE)

  # write each atom's movement vectors to a separate csv
  num_cols <- ncol(vectors)

  # Calculate the number of dataframes to create
  num_dfs <- num_cols %/% 3

  # Create a list to store the new dataframes
  new_dfs <- list()

  # Split the dataframe into the appropriate number of 3-column dataframes and store them in the list
  for (i in 1:num_dfs) {
    start_col <- (i - 1) * 3 + 1
    end_col <- start_col + 2
    new_dfs[[i]] <- vectors[,start_col:end_col]
  }
  for (i in 1:length(new_dfs)) {
    file_name <- paste0(name, '/vib_', i, '.csv')
    write.csv(new_dfs[[i]], file = file_name, row.names = FALSE)
  }
}

#' Unwrapper for feather files
#'
#' Translate stored information into human readable csv files.
#' User is guided with questions.
#' @return folder with optimized structures as xyz, and information needed for moleculaR.
#' @aliases unwRapper
#' @export
unwRapper <- function() {
  svDialogs::dlg_message('
    Please choose the output files directory, resulted from
    running extRactoR().',
    type = 'ok')$res
  output.files.dir <- svDialogs::dlg_dir()$res

  setwd(output.files.dir)

  if (dir.exists('moleculaR_csv_files')) {
    # If the folder exists, delete it
    unlink('moleculaR_csv_files', recursive = TRUE)
  }

  dir.create('moleculaR_csv_files')

  if (dir.exists('Optimized_structures_xyz')) {
    # If the folder exists, delete it
    unlink('Optimized_structures_xyz', recursive = TRUE)
  }

  dir.create('Optimized_structures_xyz')

  for (file in list.files(pattern = '.feather')) {
    name <- tools::file_path_sans_ext(file)
    moleculaR:::unwRapper.single(file)
    file.rename(name, paste0('moleculaR_csv_files',
                             '/',
                             name))
  }
  xyz_files <- list.files(list.dirs('moleculaR_csv_files/',
                                    recursive = T),
                          pattern = '.xyz',
                          full.names = T)
  names <- basename(list.dirs('moleculaR_csv_files/',
                              recursive = F))
  for (file in xyz_files) {
    invisible(file.copy(file,
              paste0('Optimized_structures_xyz/', 
                     names[[match(file, xyz_files)]], '.xyz'),
              overwrite = T
    ))
  }

  svDialogs::dlg_message(
  "
  Two new folders were created:
  
  1) Optimized_structures_xyz, which contains
  all xyz files for graphic uses and convenience.

  2) moleculaR_csv_files, which contains all of the extracted information.


  Feel free changing the xyz folder's name, but do NOT chnage the csv folder's name.",
  type = 'ok'
  )

  answer <- svDialogs::dlg_message(
    'Would you like to set the working directory to moleculaR_csv_files?

    NOTE - moleculaR() expects to be run while the working directory is
    set to "moleculaR_csv_files". it is recommended to do that now.',
    type = 'yesno')$res
  if (answer == 'yes') {
    setwd('moleculaR_csv_files/')
  }
}

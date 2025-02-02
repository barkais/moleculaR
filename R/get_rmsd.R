#' Root Mean Square Deviation calculation
#'
#' Calculate RMSD between molecular structures sharing a common substructure
#'
#' @param dir path to xyz files location
#' @param shared_indices numeric vector of atom indices defining the common substructure
#' @param reference_file optional, path to XYZ file to use as reference structure. 
#' If NULL (default), uses average of all structures
#' @return data frame with RMSD values for each structure compared to reference,
#' with structure filenames as row names  
#' @export
compute_substructure_rmsd <- function(dir = '../Optimized_structures_xyz/',
                                      shared_indices, reference_file = NULL) {
  
  original_dir <- getwd()
  on.exit(setwd(original_dir))
  setwd(dir)
  
  # Define cleanup function to remove temporary files
  cleanup_temp_files <- function() {
    temp_files <- list.files(pattern = "\\_tc.xyz$")
    if (length(temp_files) > 0) {
      unlink(temp_files)
    }
  }
  
  # Wrap the entire process in a tryCatch block
  result <- tryCatch({
    # Find all XYZ files in the current working directory
    lapply(list.files(pattern = "\\.xyz$"),
           function(x) coor.trans.file(paste(as.character(shared_indices[1:3]),
                                             collapse = ' '), x))
    xyz_files <- list.files(pattern = "\\_tc.xyz$")
    if (length(xyz_files) == 0) {
      stop("No .xyz files found in the working directory.")
    }
    
    # Read XYZ files and extract substructure coordinates
    structures <- lapply(xyz_files, function(file) {
      xyz_data <- data.table::fread(file, skip = 2, header = FALSE, col.names = c("Atom", "X", "Y", "Z"))
      xyz_sub <- xyz_data[shared_indices, ]
      as.matrix(xyz_sub[, 2:4]) # Extract coordinates (columns X, Y, Z)
    })
    
    # Validate that the selected substructure has the same size across all structures
    n_atoms <- sapply(structures, nrow)
    if (length(unique(n_atoms)) != 1) {
      stop("The selected shared substructure must have the same number of atoms across all structures.")
    }
    
    # Determine the reference structure
    if (is.null(reference_file)) {
      # Default: Compute the average structure
      reference_structure <- Reduce("+", structures) / length(structures)
    } else {
      # Use the specific reference file
      reference_file <- stringr::str_replace(reference_file, '.xyz', '_tc.xyz')
      if (!reference_file %in% xyz_files) {
        stop("The specified reference file is not found in the working directory.")
      }
      ref_index <- which(xyz_files == reference_file)
      reference_structure <- structures[[ref_index]]
    }
    
    # Calculate RMSD for each structure against the reference structure
    rmsd_values <- plyr::laply(structures, function(structure) {
      sqrt(mean(rowSums((structure - reference_structure)^2)))
    })
    
    # Create a dataframe with file names as row names
    rmsd_df <- data.frame(RMSD = rmsd_values, row.names = xyz_files)
    
    # Return the dataframe
    row.names(rmsd_df) <-  stringr::str_replace(row.names(rmsd_df), '_tc.xyz', '')
    return(rmsd_df)
  }, error = function(e) {
    # Handle errors by cleaning up temporary files and rethrowing the error
    message("An error occurred: ", e$message)
    cleanup_temp_files()
    stop(e)
  }, finally = {
    # Always clean up temporary files
    cleanup_temp_files()
  })
  
  return(result)
}



####### ----------------------------------------------------#####
####### -----------------Utility Functions------------------#####
####### ----------------------------------------------------#####

# polar.single - extracts the iso and anisotropic polarizabilty from Gaussian
#' Pull polarizability info
#' @keywords internal
#' @return A data frame with polarizabilities
polar.single <- function() {
  polar <- data.table::fread(list.files(pattern = 'polar'))
  names(polar) <- c('Iso.Polar', 'Aniso.Polar')
  return(polar)
}

# polar.df - does the same for all files in moleculaR_csv_files
#' Pull polarizability info for a set of molecules
#' @return A data frame with polarizabilities
#' @export
polar.df <- function() {
  molecules <- list.dirs(recursive = F, full.names = F)
  pol.list <- list()
  for (molecule in molecules) {
    setwd(molecule)
    pol.list[[match(molecule, molecules)]] <- polar.single()
    setwd('..')
  }
  pol.dafr <- data.frame(data.table::rbindlist(pol.list, fill = T))
  row.names(pol.dafr) <- molecules
  return(pol.dafr)
}

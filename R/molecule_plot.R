####### ----------------------------------------------------#####
####### --------------moleculaR xyz 3D Viewer---------------#####
####### ----------------------------------------------------#####

#' Generate a 3D plot of an xyz file, including atom indices
#' @param xyz_file an xyz file
#' @return plot widget (Mac and Windows only, rgl is a mess on linux)
#' @export plot_molecule
#' @export
plot_molecule <- function(xyz_file) {
  xyz <- data.table::fread(xyz_file)
  xyz$V1 <- suppressMessages(plyr::mapvalues(xyz$V1,
                                             from = c("H", "He", "Li", "Be", "B",
                                                      "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
                                                      "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
                                                      "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
                                                      "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
                                                      "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
                                                      "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
                                                      "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
                                                      "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu",
                                                      "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db",
                                                      "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
                                                      "Ts", "Og"),
                                             to = 1:118))
  atomic.numbers <- as.numeric(xyz$V1)
  
  # Plot atoms with reduced opacity
  rgl::plot3d(xyz[,2:4], col = col_and_size$V2[atomic.numbers],
              size = col_and_size$V3[atomic.numbers] * 2.2,                   # Increase atom size for visibility
              type = 's',
              alpha = 0.5,                 # Lower opacity for better text visibility
              axes = F,
              box = F)
  
  # Plot bonds with adjusted width
  bonds <- extract.connectivity(xyz_file, threshold_distance = 2.12)
  plot.edges <- unlist(lapply(1:nrow(bonds), function(i) as.list(bonds[i,])))
  rgl::segments3d(xyz[plot.edges,2:4], add = T, lwd = 2)   # Reduce bond width for clarity
  
  # Add text with improved visibility
  for (i in 1:nrow(xyz)) {
    rgl::text3d(xyz[i,2:4], texts = as.character(i),
                cex = 1.5,                  # Font size for visibility
                col = 'black',              # Text color
                fixedSize = FALSE, 
                add = T, 
                adj = c(0.5, 0.5),          # Center text on atoms
                pos = 3)                    # Position text slightly above the atoms
  }
  
  rgl::aspect3d('iso')
  rgl::rglwidget()
}


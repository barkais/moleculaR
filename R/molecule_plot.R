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
  rgl::plot3d(xyz[,2:4], col = col_and_size$V2[atomic.numbers],
              size = col_and_size$V3[atomic.numbers],
              type = 's',
              axes = F,
              box = F)
  bonds <- extract.connectivity(xyz_file, threshold_distance = 2.12)
  plot.edges <- unlist(lapply(1:nrow(bonds), function(i) as.list(bonds[i,])))
  rgl::segments3d(xyz[plot.edges,2:4], add = T, lwd = 3.5)
  rgl::text3d(xyz[,2:4], texts = 1:nrow(xyz),
              fixedSize = FALSE, add = T, adj = -1.4)
  rgl::aspect3d('iso')
  rgl::rglwidget()
}

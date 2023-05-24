####### ----------------------------------------------------#####
####### -------------moleculaR Main Functions---------------#####
####### ----------------------------------------------------#####

##### moleculaR input maker - generates Q&A, final results csv file
##### and user inputs file (.RData) for future use and documentation

#' User interface for the extraction of all possible features
#'
#' Generate user inputs for future use.
#' The input can be used directly with the moleculaR function.
#' No input is needed. User is guided with questions and answers.
#' See the full manual for rules of use and options.
#' @return csv file with all features and a input file for future
#'  use and documentation
#' @aliases moleculaR.input.maker
#' @export
moleculaR.input.maker <- function() {
  
  svDialogs::dlg_message(
    "
Welcome to moleculaR's input file maker.

Please answer and follow the instructions.",
type = 'ok')
  
  run.plot <- svDialogs::dlg_message("
It is best if you have a molecular viewer showing one of the
molecules of the set in 3D.  

Answer 'Yes' to open moleculaR's 3D viewer with a chosen molecule.

Answer 'No' if you already have an open viewer",
type = 'yesno'
  )$res
  
  if (run.plot == 'yes') {
    plot_molecule(svDialogs::dlg_open(getwd())$res)
    return(moleculaR.input.maker.rerun())
  }
  
  ### steRimol
  
  sterimol.answer <- svDialogs::dlg_message(
    "Would you like to use steRimol?",
    type = 'yesno')$res
  if (sterimol.answer == 'yes') {
    radii.system <- svDialogs::dlg_list(choices = list(
      'CPK','Pyykko'
    ), preselect = 'Pyykko',
    title = 'Please select the raddi system you wish to use.'
    )$res
    
    if (radii.system == 'CPK') CPK <- T else CPK <- F
    
    only.sub <- svDialogs::dlg_message(
      "Only account for substituent?",
      type = 'yesno'
    )$res
    
    if (only.sub == 'yes') only.sub <- T else only.sub <- F
    
    # drop is not accounted for, as it is a non standard use of steRimol.
    # See documentation and use steRimol directly to apply it.
    
    input.vector <- svDialogs::dlg_input(
      "Please enter atom pairs (min of one pair).
Separate atom pair by a single comma.

If you choose to enter more than one pair,
separate your answers with a space.

example: to use atom pairs 1,2 and 3,4
as primary axes, our answer should be - 1,2 3,4,

Note: Directionality matters. 1,2 and 2,1 are different."
    )$res
    
    input.vector <- stringr::str_split(input.vector, ' ')
    input.vector <- lapply(input.vector,
                           function(x) gsub(',', ' ', x))[[1]]
    
    sterimol.inputs <- list(input.vector,
                            CPK,
                            only.sub)
    names(sterimol.inputs) <- c('input.vector',
                                'CPK',
                                'only.sub')
  } else {
    sterimol.inputs <- NA
  }
  
  ### NBO
  
  nbo.answer <- svDialogs::dlg_message(
    "Would you like to use NBO charges?",
    type = 'yesno')$res
  if (nbo.answer == 'yes') {
    
    atom_indices <- svDialogs::dlg_input(
      "Please enter atoms, separated by a comma"
    )$res
    
    atom_indices <- gsub(',', ' ', atom_indices)
    
    diff.answer <- svDialogs::dlg_message(
      "Would you like to use NBO charge differences?",
      type = 'yesno')$res
    
    if (diff.answer == 'yes') {
      
      difference_indices <- svDialogs::dlg_input(
        "Please enter atom pairs, atoms separated by a comma
and pairs by a space"
      )$res
      
      difference_indices <- gsub(',', ' ', difference_indices)
    } else {
      difference_indices <- '-'
    }
    
    nbo.inputs <- list(atom_indices,
                       difference_indices)
    names(nbo.inputs) <- c('atom_indices',
                           'difference_indices')
  } else {
    nbo.inputs <- NA
  }
  
  
  ### Dipole
  
  dipole.answer <- svDialogs::dlg_message(
    "Would you like to use dipole moment(s)?",
    type = 'yesno')$res
  
  if (dipole.answer == 'yes') {
    change.coor.sys <- svDialogs::dlg_message(
      "Would you like to perform any change of coordinates system?

This is a must if you wish to use any dipole manipulation.",
type = 'yesno'
    )$res
    
    if (change.coor.sys == 'yes') {
      coor_atoms <- svDialogs::dlg_input(
        "
  Please enter three atom indices.
  The first being the origin atom,
  the y direction atom, and the xy plane atom.
  The answer should be 3 numbers, separated by a comma (e.g. 1,2,3).
  
  If you choose to enter more than one triad,
  separate your answers with a space.
        "
      )$res
      
      coor_atoms <- stringr::str_split(coor_atoms, ' ')
      coor_atoms <- lapply(coor_atoms,
                           function(x) gsub(',', ' ', x))[[1]]
      
      center_of_substructure <- svDialogs::dlg_message(
        "Would you like to use the centeroid
of a substructure as the origin of the dipole components?",
type = 'yesno'
      )$res
      
      if (center_of_substructure == 'no') {
        center_of_substructure <- F
        subunits_input_vector <- NULL
        center_of_mass <- F
      } else {
        center_of_mass <- F
        center_of_substructure <- T
        subunits_input_vector <- svDialogs::dlg_input(
          "Please enter the subunit's atoms (min of one subunit).
Separate atom pair by a single comma.

If you choose to enter more than one subunit,
separate your answers with a space.

example: if we wish to use atom pairs 1,2,3,4 and 5,6,7,8
as primary axes, our answer should be - 1,2,3,4 and 5,6,7,8"
        )$res
        
        subunits_input_vector <- stringr::str_split(subunits_input_vector, ' ')
        subunits_input_vector <- lapply(subunits_input_vector,
                                        function(x) gsub(',', ' ', x))[[1]]
      }
    } else {
      coor_atoms <- ''
      center_of_mass <- F
      center_of_substructure <- F
      subunits_input_vector <- NULL
    }
    
    dipole.inputs <- list(coor_atoms,
                          center_of_substructure,
                          subunits_input_vector,
                          center_of_mass)
    names(dipole.inputs) <- c('coor_atoms',
                              'center_of_substructure',
                              'subunits_input_vector',
                              'center_of_mass')
  } else {
    dipole.inputs <- NA
  }
  
  
  ###### Vibrations
  
  ### Bond Vibrations
  
  bond.vibrations.answer <- svDialogs::dlg_message(
    "Would you like to use bond vibrations?",
    type = 'yesno')$res
  
  if (bond.vibrations.answer == 'yes') {
    atom_pairs <- svDialogs::dlg_input(
      "
Please enter atom pairs. 

Each pair should be separated by a comma
and different pairs should be separated by a space (e.g. 1,2 3,4). 

Note: there must be a bond between the requested atoms.
"
    )$res
    atom_pairs <- stringr::str_split(atom_pairs, ' ')
    atom_pairs <- paste(lapply(atom_pairs,
                               function(x) gsub(',', ' ', x))[[1]],
                        collapse =  ' ')
    stretch.inputs <- list(atom_pairs)
    names(stretch.inputs) <- 'atom_pairs'
  } else {
    stretch.inputs <- NA
  }
  
  ### Ring Vibrations
  
  ring.vibrations.answer <- svDialogs::dlg_message(
    "Would you like to use ring vibrations?",
    type = 'yesno')$res
  
  if (ring.vibrations.answer == 'yes') {
    inputs_vector <- svDialogs::dlg_input(
      "

Please enter sets of ring atoms by order (para, primary, o, o, m, m).

You may input more than one set. All atoms in a set are to be separated
by a comma, and different sets by a space.


")$res
    
    inputs_vector <- stringr::str_split(inputs_vector, ' ')
    inputs_vector <- lapply(inputs_vector,
                            function(x) gsub(',', ' ', x))[[1]]
    ring.inputs <- list(inputs_vector)
    names(ring.inputs) <- 'inputs_vector'
  } else {
    ring.inputs <- NA
  }
  
  ### Bending Vibrations
  
  bending.vibrations.answer <- svDialogs::dlg_message(
    "Would you like to use bending vibrations?",
    type = 'yesno')$res
  
  if (bending.vibrations.answer == 'yes') {
    atom_pairs <- svDialogs::dlg_input(
      "
Please enter atom pairs. 

The pair should be pf the atoms that move, and not their center atom. 

Each pair should be separated by a comma
and different pairs should be separated by a space (e.g. 1,2 3,4)


"
    )$res
    atom_pairs <- stringr::str_split(atom_pairs, ' ')
    atom_pairs <- paste(lapply(atom_pairs,
                               function(x) gsub(',', ' ', x))[[1]],
                        collapse =  ' ')
    bend.inputs <- list(atom_pairs)
    names(bend.inputs) <- 'atom_pairs'
  } else {
    bend.inputs <- NA
  }
  
  ### Angles
  
  angles.answer <- svDialogs::dlg_message(
    "Would you like to use angles and/or dihedrals?",
    type = 'yesno')$res
  
  if (angles.answer == 'yes') {
    atom_sets <- svDialogs::dlg_input(
      "
Please enter atom triads for angles or quartets for dihedrals.
Each atom in the triad/quartet should be separated by a comma,
and different sets should be separated by a space (e.g. 1,2,3 1,2,3,4).


"
    )$res
    
    inputs_vector <- stringr::str_split(atom_sets, ' ')
    inputs_vector <- lapply(inputs_vector,
                            function(x) gsub(',', ' ', x))[[1]]
    mol.angle.inputs <- list(inputs_vector)
    names(mol.angle.inputs) <- 'inputs_vector'
  } else {
    mol.angle.inputs <- NA
  }
  
  ### Distances
  
  distance.answer <- svDialogs::dlg_message(
    "Would you like to use atom distances?",
    type = 'yesno')$res
  
  if (distance.answer == 'yes') {
    atom_sets <- svDialogs::dlg_input(
      "
Please enter atom pairs for which you wish to the distance between.
Each atom should be separated by a comma,
and different pairs should be separated by a space (e.g. 1,2,3 1,2,3,4).


"
    )$res
    
    inputs_vector <- stringr::str_split(atom_sets, ' ')
    inputs_vector <- paste(lapply(inputs_vector,
                                  function(x) gsub(',', ' ', x))[[1]],
                           collapse =  ' ')
    distances.inputs <- list(inputs_vector)
    names(distances.inputs) <- 'inputs_vector'
  } else {
    distances.inputs <- NA
  }
  
  ### Polarizability
  
  ### Polarizability
  
  polariz.answer <- svDialogs::dlg_message(
    "Would you like to use iso an anisotropic polarizabilty?",
    type = 'yesno')$res
  
  polar.inputs <- list(polariz.answer)
  names(polar.inputs) <- 'polariz.answer'
  
  
  
  #####################################################################
  ###################   Output File Preparation   #####################
  #####################################################################
  
  inputs.list <- list(sterimol.inputs,
                      nbo.inputs,
                      dipole.inputs,
                      stretch.inputs,
                      ring.inputs,
                      bend.inputs,
                      mol.angle.inputs,
                      distances.inputs,
                      polar.inputs)
  names(inputs.list) <- c('steRimol',
                          'NBO',
                          'Dipole',
                          'Bond Vibs',
                          'Ring Vibs',
                          'Bend Vibs',
                          'Angles',
                          'Distances',
                          'Polarizability')
  
  name.of.file <- svDialogs::dlg_input(
    "
How do you want to call the input file?

Do not add a file extension."
  )$res
  
  svDialogs::dlg_message(
    "
Please choose where to save the input file.",
type = 'ok'
  )
  where.to.save <- svDialogs::dlg_dir()$res
  saveRDS(inputs.list,paste0(where.to.save,
                             '/',
                             name.of.file,
                             '.RData'))
}

#' User interface for the extraction of all possible features - looper
#'
#' Generate user inputs for future use.
#' The input can be used directly with the moleculaR function.
#' No input is needed. User is guided with questions and answers.
#' See the full manual for rules of use and options.
#' @return csv file with all features and a input file for future
#'  use and documentation
#' @aliases moleculaR.input.maker.rerun
#' @keywords internal

moleculaR.input.maker.rerun <- function() {
  
  ### steRimol
  
  sterimol.answer <- svDialogs::dlg_message(
    "Would you like to use steRimol?",
    type = 'yesno')$res
  if (sterimol.answer == 'yes') {
    radii.system <- svDialogs::dlg_list(choices = list(
      'CPK','Pyykko'
    ), preselect = 'Pyykko',
    title = 'Please select the raddi system you wish to use.'
    )$res
    
    if (radii.system == 'CPK') CPK <- T else CPK <- F
    
    only.sub <- svDialogs::dlg_message(
      "Only account for substituent?",
      type = 'yesno'
    )$res
    
    if (only.sub == 'yes') only.sub <- T else only.sub <- F
    
    # drop is not accounted for, as it is a non standard use of steRimol.
    # See documentation and use steRimol directly to apply it.
    
    input.vector <- svDialogs::dlg_input(
      "Please enter atom pairs (min of one pair).
Separate atom pair by a single comma.

If you choose to enter more than one pair,
separate your answers with a space.

example: to use atom pairs 1,2 and 3,4
as primary axes, our answer should be - 1,2 3,4,

Note: Directionality matters. 1,2 and 2,1 are different."
    )$res
    
    input.vector <- stringr::str_split(input.vector, ' ')
    input.vector <- lapply(input.vector,
                           function(x) gsub(',', ' ', x))[[1]]
    
    sterimol.inputs <- list(input.vector,
                            CPK,
                            only.sub)
    names(sterimol.inputs) <- c('input.vector',
                                'CPK',
                                'only.sub')
  } else {
    sterimol.inputs <- NA
  }
  
  ### NBO
  
  nbo.answer <- svDialogs::dlg_message(
    "Would you like to use NBO charges?",
    type = 'yesno')$res
  if (nbo.answer == 'yes') {
    
    atom_indices <- svDialogs::dlg_input(
      "Please enter atoms, separated by a comma"
    )$res
    
    atom_indices <- gsub(',', ' ', atom_indices)
    
    diff.answer <- svDialogs::dlg_message(
      "Would you like to use NBO charge differences?",
      type = 'yesno')$res
    
    if (diff.answer == 'yes') {
      
      difference_indices <- svDialogs::dlg_input(
        "Please enter atom pairs, atoms separated by a comma
and pairs by a space"
      )$res
      
      difference_indices <- gsub(',', ' ', difference_indices)
    } else {
      difference_indices <- '-'
    }
    
    nbo.inputs <- list(atom_indices,
                       difference_indices)
    names(nbo.inputs) <- c('atom_indices',
                           'difference_indices')
  } else {
    nbo.inputs <- NA
  }
  
  
  ### Dipole
  
  dipole.answer <- svDialogs::dlg_message(
    "Would you like to use dipole moment(s)?",
    type = 'yesno')$res
  
  if (dipole.answer == 'yes') {
    change.coor.sys <- svDialogs::dlg_message(
      "Would you like to perform any change of coordinates system?

This is a must if you wish to use any dipole manipulation.",
type = 'yesno'
    )$res
    
    if (change.coor.sys == 'yes') {
      coor_atoms <- svDialogs::dlg_input(
        "
  Please enter three atom indices.
  The first being the origin atom,
  the y direction atom, and the xy plane atom.
  The answer should be 3 numbers, separated by a comma (e.g. 1,2,3).
  
  If you choose to enter more than one triad,
  separate your answers with a space.
        "
      )$res
      
      coor_atoms <- stringr::str_split(coor_atoms, ' ')
      coor_atoms <- lapply(coor_atoms,
                           function(x) gsub(',', ' ', x))[[1]]
      
      center_of_substructure <- svDialogs::dlg_message(
        "Would you like to use the centeroid
of a substructure as the origin of the dipole components?",
type = 'yesno'
      )$res
      
      if (center_of_substructure == 'no') {
        center_of_substructure <- F
        subunits_input_vector <- NULL
        center_of_mass <- F
      } else {
        center_of_mass <- F
        center_of_substructure <- T
        subunits_input_vector <- svDialogs::dlg_input(
          "Please enter the subunit's atoms (min of one subunit).
Separate atom pair by a single comma.

If you choose to enter more than one subunit,
separate your answers with a space.

example: if we wish to use atom pairs 1,2,3,4 and 5,6,7,8
as primary axes, our answer should be - 1,2,3,4 and 5,6,7,8"
        )$res
        
        subunits_input_vector <- stringr::str_split(subunits_input_vector, ' ')
        subunits_input_vector <- lapply(subunits_input_vector,
                                        function(x) gsub(',', ' ', x))[[1]]
      }
    } else {
      coor_atoms <- ''
      center_of_mass <- F
      center_of_substructure <- F
      subunits_input_vector <- NULL
    }
    
    dipole.inputs <- list(coor_atoms,
                          center_of_substructure,
                          subunits_input_vector,
                          center_of_mass)
    names(dipole.inputs) <- c('coor_atoms',
                              'center_of_substructure',
                              'subunits_input_vector',
                              'center_of_mass')
  } else {
    dipole.inputs <- NA
  }
  
  
  ###### Vibrations
  
  ### Bond Vibrations
  
  bond.vibrations.answer <- svDialogs::dlg_message(
    "Would you like to use bond vibrations?",
    type = 'yesno')$res
  
  if (bond.vibrations.answer == 'yes') {
    atom_pairs <- svDialogs::dlg_input(
      "
Please enter atom pairs. Each pair should be separated by a comma
and different pairs should be separated by a space (e.g. 1,2 3,4)

"
    )$res
    atom_pairs <- stringr::str_split(atom_pairs, ' ')
    atom_pairs <- paste(lapply(atom_pairs,
                               function(x) gsub(',', ' ', x))[[1]],
                        collapse =  ' ')
    stretch.inputs <- list(atom_pairs)
    names(stretch.inputs) <- 'atom_pairs'
  } else {
    stretch.inputs <- NA
  }
  
  ### Ring Vibrations
  
  ring.vibrations.answer <- svDialogs::dlg_message(
    "Would you like to use ring vibrations?",
    type = 'yesno')$res
  
  if (ring.vibrations.answer == 'yes') {
    inputs_vector <- svDialogs::dlg_input(
      "

Please enter sets of ring atoms by order (para, primary, o, o, m, m).

You may input more than one set. All atoms in a set are to be separated
by a comma, and different sets by a space.


")$res
    
    inputs_vector <- stringr::str_split(inputs_vector, ' ')
    inputs_vector <- lapply(inputs_vector,
                            function(x) gsub(',', ' ', x))[[1]]
    ring.inputs <- list(inputs_vector)
    names(ring.inputs) <- 'inputs_vector'
  } else {
    ring.inputs <- NA
  }
  
  ### Bending Vibrations
  
  bending.vibrations.answer <- svDialogs::dlg_message(
    "Would you like to use bending vibrations?",
    type = 'yesno')$res
  
  if (bending.vibrations.answer == 'yes') {
    atom_pairs <- svDialogs::dlg_input(
      "
Please enter atom pairs. Each pair should be separated by a comma
and different pairs should be separated by a space (e.g. 1,2 3,4)


"
    )$res
    atom_pairs <- stringr::str_split(atom_pairs, ' ')
    atom_pairs <- paste(lapply(atom_pairs,
                               function(x) gsub(',', ' ', x))[[1]],
                        collapse =  ' ')
    bend.inputs <- list(atom_pairs)
    names(bend.inputs) <- 'atom_pairs'
  } else {
    bend.inputs <- NA
  }
  
  ### Angles
  
  angles.answer <- svDialogs::dlg_message(
    "Would you like to use angles and/or dihedrals?",
    type = 'yesno')$res
  
  if (angles.answer == 'yes') {
    atom_sets <- svDialogs::dlg_input(
      "
Please enter atom triads for angles or quartets for dihedrals.
Each atom in the triad/quartet should be separated by a comma,
and different sets should be separated by a space (e.g. 1,2,3 1,2,3,4).


"
    )$res
    
    inputs_vector <- stringr::str_split(atom_sets, ' ')
    inputs_vector <- lapply(inputs_vector,
                            function(x) gsub(',', ' ', x))[[1]]
    mol.angle.inputs <- list(inputs_vector)
    names(mol.angle.inputs) <- 'inputs_vector'
  } else {
    mol.angle.inputs <- NA
  }
  
  ### Distances
  
  distance.answer <- svDialogs::dlg_message(
    "Would you like to use atom distances?",
    type = 'yesno')$res
  
  if (distance.answer == 'yes') {
    atom_sets <- svDialogs::dlg_input(
      "
Please enter atom pairs for which you wish to the distance between.
Each atom should be separated by a comma,
and different pairs should be separated by a space (e.g. 1,2,3 1,2,3,4).


"
    )$res
    
    inputs_vector <- stringr::str_split(atom_sets, ' ')
    inputs_vector <- paste(lapply(inputs_vector,
                                  function(x) gsub(',', ' ', x))[[1]],
                           collapse =  ' ')
    distances.inputs <- list(inputs_vector)
    names(distances.inputs) <- 'inputs_vector'
  } else {
    distances.inputs <- NA
  }
  
  ### Polarizability
  
  polariz.answer <- svDialogs::dlg_message(
    "Would you like to use iso an anisotropic polarizabilty?",
    type = 'yesno')$res
  
  polar.inputs <- list(polariz.answer)
  names(polar.inputs) <- 'polariz.answer'
  
  
  
  
  #####################################################################
  ###################   Output File Preparation   #####################
  #####################################################################
  
  inputs.list <- list(sterimol.inputs,
                      nbo.inputs,
                      dipole.inputs,
                      stretch.inputs,
                      ring.inputs,
                      bend.inputs,
                      mol.angle.inputs,
                      distances.inputs,
                      polar.inputs)
  names(inputs.list) <- c('steRimol',
                          'NBO',
                          'Dipole',
                          'Bond Vibs',
                          'Ring Vibs',
                          'Bend Vibs',
                          'Angles',
                          'Distances',
                          'Polarizability')
  
  name.of.file <- svDialogs::dlg_input(
    "
How do you want to call the input file?

Do not add a file extension."
  )$res
  
  svDialogs::dlg_message(
    "
Please choose where to save the input file.",
type = 'ok'
  )
  where.to.save <- svDialogs::dlg_dir()$res
  saveRDS(inputs.list,paste0(where.to.save,
                             '/',
                             name.of.file,
                             '.RData'))
} 
##### moleculaR - generates final results csv file
##### from user inputs file (.RData), after running GUI for the first time

#' User function for the extraction of all possible features
#' based on a ready input file
#'
#' Calculate all features by applying to an input file, saving the final
#'  data set in a chosen location.
#' @param input_file moleculaR's input file
#' @return csv file with all features and a input file for future
#'  use and documentation
#' @aliases moleculaR.input
#' @export
moleculaR <- function(input_file = NULL) {

  home <- getwd()
  #####################################################################
  #############   Feature Computation and Extraction   ################
  #####################################################################

  if (is.null(input_file)) {
    input_file <- readRDS(svDialogs::dlg_open(title =
                                              "Please choose an .RData input file",
                                            getwd())$res)
  } else {
    input_file <- readRDS(input_file)
  }
  results <- list()

  ### steRimol

  if (!any(is.na(input_file$steRimol))) {
    sterimol.result <- list(steRimol.multi(input_file$steRimol$input.vector,
                                               input_file$steRimol$CPK,
                                               input_file$steRimol$only.sub))
    results <- c(results, sterimol.result)
  }

  ### NBO

  if (!any(is.na(input_file$NBO))) {
    nbo.result <- list(nbo.df(input_file$NBO$atom_indices,
                                   input_file$NBO$difference_indices))
    results <- c(results, nbo.result)
  }

  ### Dipole

  if (!any(is.na(input_file$Dipole))) {
    dipole.result <- list(dip.gaussian.multi(input_file$Dipole$coor_atoms,
                                             input_file$Dipole$center_of_substructure,
                                             input_file$Dipole$subunits_input_vector,
                                             input_file$Dipole$center_of_mass))
    results <- c(results, dipole.result)
  }

  ###### Vibrations

  ### Bond Vibrations

  if (!is.na(input_file$`Bond Vibs`)) {
    bond.vibrations.result <- list(stretch.vib.df(input_file$`Bond Vibs`$atom_pairs))
    results <- c(results, bond.vibrations.result)
  }

  ### Ring Vibrations

  if (!is.na(input_file$`Ring Vibs`)) {
    ring.vibrations.result <- list(ring.vib.multi(input_file$`Ring Vibs`$inputs_vector))
    results <- c(results, ring.vibrations.result)
  }

  ### Bending Vibrations

  if (!is.na(input_file$`Bend Vibs`)) {
    bend.vibrations.result <- list(bend.vib.df(input_file$`Bend Vibs`$atom_pairs))
    results <- c(results, bend.vibrations.result)
  }

  ### Angles

  if (!is.na(input_file$Angles)) {
    angles.result <- list(mol.angles.multi(input_file$Angles$inputs_vector))
    results <- c(results, angles.result)
  }

  ### Distances

  if (!is.na(input_file$Distances)) {
    distances.result <- list(atoms.distance.df(input_file$Distances$inputs_vector))
    results <- c(results, distances.result)
  }

  ### Polarizability

  if (input_file$Polarizability[1] == 'yes') {
    polar.result <- list(polar.df())
    results <- c(results, polar.result)
  } 



  #####################################################################
  ###################   Output File Preparation   #####################
  #####################################################################



  results.df <- do.call(cbind, results)
  name.of.csv <- svDialogs::dlg_input(
    "
How do you want to call the features csv file?

Do not add a file extension."
  )$res

  svDialogs::dlg_message(
    "
Please choose where to save the csv file, and inputs.",
type = 'ok'
  )
  where.to.save <- svDialogs::dlg_dir()$res
  write.csv(results.df,
            paste0(where.to.save,
                   '/',
                   name.of.csv,
                   '.csv'))

  on.exit(setwd(home))

}


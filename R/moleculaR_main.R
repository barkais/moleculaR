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
#' @export
#' @import shinythemes
#' @import shiny
moleculaR.Input.Maker <- function() {
  # UI
  ui <- fluidPage(
    theme = shinythemes::shinytheme("lumen"),  # Apply a modern theme
    
    tags$head(
      tags$script(HTML("
        $(document).on('click', '[id$=InfoBtn]', function() {
          var infoBoxId = $(this).attr('id').replace('Btn', 'Box');
          $('#' + infoBoxId).toggle();
        });
      ")),
      tags$style(HTML("
        .info-box {
          display: none;
          background-color: #f9f9f9;
          border-left: 4px solid #5bc0de;
          padding: 10px;
          margin-top: 10px;
        }
        .section-title {
          margin-top: 20px;
          font-weight: bold;
          color: #337ab7;
        }
        .section-container {
          padding: 15px;
          border: 1px solid #ddd;
          border-radius: 5px;
          background-color: #ffffff;
          margin-bottom: 20px;
        }
        .info-btn {
          font-size: 14px;
          color: #5bc0de;
          margin-left: 5px;
          cursor: pointer;
        }
        .section-divider {
          margin-top: 30px;
          margin-bottom: 30px;
          height: 2px;
          background-color: #eee;
        }
      "))
    ),
    
    titlePanel("moleculaR - Input File Maker"),
    
    fluidRow(
      column(4, 
             div(class = "section-container",
                 h3("Sterimol", span("\u2139\ufe0f", id = "sterimolInfoBtn", class = "info-btn")),
                 div(id = "sterimolInfoBox", class = "info-box", style = "white-space: pre-wrap;",
                     "Enter pairs of integers separated by a space,
  '_', tab, ':', or ';'
  For example:
  '1,2 3,4' or '1,2_3,4' or '1,2;3,4'"),
  checkboxInput("sterimol", "Sterimol", value = FALSE),
  conditionalPanel(
    condition = "input.sterimol == true",
    textInput("axes", "Axis (or axes):", ""),
    selectInput("radii", "Radii:", choices = c("CPK", "Pyykko")),
    checkboxInput("onlySubstituent", "Only substituent?", value = TRUE),
    checkboxInput("dropAtoms", "Drop any atoms?", value = FALSE),
    conditionalPanel(
      condition = "input.dropAtoms == true",
      textInput("atoms", "Enter atoms to drop:", "")
    )
  )
             )
      ),
  
  column(4, 
         div(class = "section-container",
             h3("Dipole Moment", span("\u2139\ufe0f", id = "dipoleInfoBtn", class = "info-btn")),
             div(id = "dipoleInfoBox", class = "info-box", style = "white-space: pre-wrap;",
                 "Enter vectors of integers separated by a space,
  '_', tab, ':', or ';'
  For example:
  '1,2,3,4 5,6,7' or '1,2,3,4_5,6,7'"),
  checkboxInput("dipoleMoment", "Dipole Moment", value = FALSE),
  conditionalPanel(
    condition = "input.dipoleMoment == true",
    textInput("dipoleInput", "Enter vectors of atom indices:", "")
  )
         )
  ),
  
  column(4, 
         div(class = "section-container",
             h3("Charges", span("\u2139\ufe0f", id = "chargesInfoBtn", class = "info-btn")),
             div(id = "chargesInfoBox", class = "info-box", style = "white-space: pre-wrap;",
                 "Atom indices input line:
  Enter a vector of integers. Separate indices with commas.
  For example:
  '1,2,3,4'
  
  Differences input line: 
  Enter pairs of integers that all appear in the first input line. 
  Pairs should be separated by space, '_', tab, ':', or ';'. 
  For example:
  '1,2 3,4' or '1,2_3,4' or '1,2;3,4'"),
  checkboxInput("charges", "Charges", value = FALSE),
  conditionalPanel(
    condition = "input.charges == true",
    textInput("atomIndices", "Atom indices:", ""),
    textInput("differences", "Differences:", ""),
    checkboxGroupInput("chargeMethods", "Select charge methods:", 
                       choices = list("NPA" = "NPA", "Hirshfeld" = "Hirshfeld", "CM5" = "CM5"))
  )
         )
  )
    ),
  
  div(class = "section-divider"),  # Divider for clarity
  
  fluidRow(
    column(4, 
           div(class = "section-container",
               h3("Vibrations", span("\u2139\ufe0f", id = "vibrationsInfoBtn", class = "info-btn")),
               div(id = "vibrationsInfoBox", class = "info-box", style = "white-space: pre-wrap;",
                   "Bond vibrations input line: 
  Enter pairs of integers separated by space, '_', tab, ':', or ';'. 
  For example:
  '1,2 3,4' or '1,2_3,4' or '1,2;3,4'.
  
  Ring vibrations input line: 
  Enter integers separated by space, '_', tab, ':', or ';'. 
  For example:
  '1 2 3' or '1_2_3' or '1;2;3'.
  
  Bending vibrations input line: 
  Enter pairs of integers separated by space, '_', tab, ':', or ';'. 
  For example:
  '1,2 3,4' or '1,2_3,4' or '1,2;3,4'."),
  checkboxInput("vibrations", "Vibrations", value = FALSE),
  conditionalPanel(
    condition = "input.vibrations == true",
    checkboxInput("bondVibrations", "Bond Vibrations", value = FALSE),
    conditionalPanel(
      condition = "input.bondVibrations == true",
      textInput("bondVibInput", "Enter bonded atom pairs:", "")
    ),
    checkboxInput("ringVibrations", "Ring Vibrations", value = FALSE),
    conditionalPanel(
      condition = "input.ringVibrations == true",
      textInput("ringVibInput", "Enter the primary atom(s) index:", "")
    ),
    checkboxInput("bendVibrations", "Bend Vibrations", value = FALSE),
    conditionalPanel(
      condition = "input.bendVibrations == true",
      textInput("bendVibInput", "Enter pairs of atoms (terminal, moving atoms):", "")
    )
  )
           )
    ),
  
  column(4, 
         div(class = "section-container",
             h3("Geometric Measurements", span("\u2139\ufe0f", id = "additionalInfoBtn", class = "info-btn")),
             div(id = "additionalInfoBox", class = "info-box", style = "white-space: pre-wrap;",
                 "Angles input line: 
  Enter vectors of 3 or 4 integers separated by space, '_', tab, ':', or ';'. 
  For example:
  '1,2,3 4,5,6,7' or '1,2,3_4,5,6,7' or '1,2,3;4,5,6,7'.
  
  Distances input line: 
  Enter pairs of integers separated by space, '_', tab, ':', or ';'. 
  For example:
  '1,2 3,4' or '1,2_3,4' or '1,2;3,4'."),
  checkboxInput("angles", "Angles (dihedral and between bonds)", value = FALSE),
  conditionalPanel(
    condition = "input.angles == true",
    textInput("anglesInput", "Enter vectors of indices(3 for angle and 4 for dihedrals):", "")
  ),
  checkboxInput("distances", "Distances (and bond lengths)", value = FALSE),
  conditionalPanel(
    condition = "input.distances == true",
    textInput("distancesInput", "Enter pairs of atoms:", "")
  )
         )
  )
  ),
  
  div(class = "section-divider"),  # Divider for clarity
  
  fluidRow(
    column(12, 
           div(class = "section-container",
               h3("Molecule Visualization"),
               checkboxInput("molViz", "Molecule Visualization", value = FALSE),
               conditionalPanel(
                 condition = "input.molViz == true",
                 fileInput("xyzFile", "Choose an XYZ file", accept = c(".xyz")),
                 rglwidgetOutput("molPlot")
               )
           )
    )
  ),
  
  fluidRow(
    column(12,
           div(class = "section-container",
               downloadButton("downloadFile", "Download input file"),
               textOutput("outputText")
           )
    )
  )
  )
  
  # Server
  server <- function(input, output) {
    output$outputText <- renderText({
      # Initialize an empty output string
      output_data <- ""
      
      if (input$sterimol) {
        # Gather Sterimol input values
        axes <- input$axes
        radii <- input$radii
        only_substituent <- ifelse(input$onlySubstituent, "Yes", "No")
        drop_atoms <- ifelse(input$dropAtoms, input$atoms, "No")
        
        # Append Sterimol data to output
        output_data <- paste0(
          output_data,
          "Sterimol:\n",
          "Axis (or axes): ", axes, "\n",
          "Radii: ", radii, "\n",
          "Only substituent?: ", only_substituent, "\n",
          "Drop any atoms?: ", drop_atoms, "\n"
        )
      }
      
      if (input$dipoleMoment) {
        # Gather Dipole Moment input value
        dipole_value <- input$dipoleInput
        
        # Append Dipole Moment data to output
        output_data <- paste0(
          output_data,
          "Dipole Moment: ", dipole_value, "\n"
        )
      }
      
      if (input$charges) {
        # Gather Charges input values
        atom_indices <- input$atomIndices
        differences <- input$differences
        charge_methods <- paste(input$chargeMethods, collapse = ", ")
        
        # Append Charges data to output
        output_data <- paste0(
          output_data,
          "Charges:\n",
          "Atom indices: ", atom_indices, "\n",
          "Differences: ", differences, "\n",
          "Charge Methods: ", charge_methods, "\n"
        )
      }
      
      if (input$vibrations) {
        # Gather Vibrations input values
        bond_vib <- ifelse(input$bondVibrations, input$bondVibInput, "No")
        ring_vib <- ifelse(input$ringVibrations, input$ringVibInput, "No")
        bend_vib <- ifelse(input$bendVibrations, input$bendVibInput, "No")
        
        # Append Vibrations data to output
        output_data <- paste0(
          output_data,
          "Vibrations:\n",
          "Bond Vibrations: ", bond_vib, "\n",
          "Ring Vibrations: ", ring_vib, "\n",
          "Bend Vibrations: ", bend_vib, "\n"
        )
      }
      
      if (input$angles) {
        angles_data <- input$anglesInput
        output_data <- paste0(output_data, "Angles: ", angles_data, "\n")
      }
      
      if (input$distances) {
        distances_data <- input$distancesInput
        output_data <- paste0(output_data, "Distances: ", distances_data, "\n")
      }
      
      # Save the output to a file
      writeLines(output_data, "output.txt")
      
      # Return output data to display in the app (optional)
      return(output_data)
    })
    
    output$downloadFile <- downloadHandler(
      filename = function() {
        "moleculaR_inputs.txt"
      },
      content = function(file) {
        file.copy("output.txt", file)
        
        # Delete the "output.txt" file after copying
        file.remove("output.txt")
      }
    )
    
    observeEvent(input$xyzFile, {
      req(input$xyzFile)
      
      # Use the plot_molecule function to plot the molecule from the selected file
      plot_molecule(input$xyzFile$datapath)
      
      # Convert the RGL plot to a widget for Shiny
      output$molPlot <- renderRglwidget({
        rglwidget()
      })
    })
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
}

##### Input.file.Parser - translates text form input files generated with
##### moleculaR.Input.Maker. 

#' Text Input file parser to moleculaR's main function
#'
#' Transform into list form.
#' @keywords internal
#' @return list form of input file
parse_txt_inputfile_to_list <- function(input_file_path) {
  # Initialize the overall list
  Input.File <- list()
  
  # Read in the text file
  lines <- readLines(input_file_path)
  
  # Variables to keep track of the current primary list and secondary list
  current_primary <- NULL
  
  # Helper function to convert "Yes"/"No" to TRUE/FALSE
  convert_to_logical <- function(value) {
    if (tolower(value) == "yes") {
      return(TRUE)
    } else if (tolower(value) == "no") {
      return(FALSE)
    } else {
      return(value)
    }
  }
  
  # Helper function to split characters and replace commas with spaces
  process_characters <- function(value) {
    # Split at spaces, tabs, underscores, colons, or semicolons
    value <- strsplit(value, "[ \t_:;]+")[[1]]
    # Remove leading and trailing spaces and replace commas with spaces
    value <- trimws(gsub(",", " ", value))
    return(value)
  }
  
  # Loop through each line
  for (line in lines) {
    # Trim any leading or trailing whitespace
    line <- trimws(line)
    
    # Skip empty lines
    if (nchar(line) == 0) next
    
    # Check if the line ends with a colon (indicates a primary list)
    if (grepl(":$", line)) {
      # Extract the primary list name (remove the colon)
      primary_name <- sub(":", "", line)
      # Initialize a new primary list in Input.File
      Input.File[[primary_name]] <- list()
      current_primary <- primary_name
    } else if (!is.null(current_primary)) {
      # Split the line at the first colon
      parts <- strsplit(line, ":", fixed = TRUE)[[1]]
      if (length(parts) == 2) {
        # Extract the secondary list name and value
        secondary_name <- trimws(parts[1])
        value <- convert_to_logical(trimws(parts[2]))
        # Process characters if the value is not logical (TRUE/FALSE)
        if (!is.logical(value)) {
          value <- process_characters(value)
        }
        
        # Handle special cases
        if (current_primary == "Sterimol" && secondary_name == "Radii" && value == "CPK") {
          value <- TRUE
        }
        
        # Store the value in the corresponding primary list
        Input.File[[current_primary]][[secondary_name]] <- value
      }
    }
  }
  
  # Ensure "Dipole Moment", "Angles", and "Distances" are treated as primary lists
  for (special_list in c("Dipole Moment", "Angles", "Distances")) {
    if (!is.null(Input.File[[special_list]])) {
      next  # Already handled
    }
    
    # Look for the special primary list in the remaining text lines
    for (line in lines) {
      if (startsWith(line, special_list)) {
        # Extract the value after the colon
        value <- sub(paste0(special_list, ":"), "", line)
        value <- convert_to_logical(trimws(value))
        if (!is.logical(value)) {
          value <- process_characters(value)
        }
        Input.File[[special_list]] <- value
      }
    }
  }
  
  # Remove specific nested elements
  Input.File$Sterimol$`Dipole Moment` <- NULL
  Input.File$Vibrations$Angles <- NULL
  Input.File$Vibrations$Distances <- NULL
  
  return(Input.File)
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
    input_file <- parse_txt_inputfile_to_list(svDialogs::dlg_open(title =
                                              "Please choose a moleculaR.Input.Maker generated input file",
                                            getwd())$res)
  } else {
    input_file <- parse_txt_inputfile_to_list(input_file)
  }
  results <- list()

  ### steRimol

  if ('Sterimol' %in% names(input_file)) {
    sterimol.result <- list(steRimol.multi(input_file$Sterimol$`Axis (or axes)`,
                                           input_file$Sterimol$Radii,
                                           input_file$Sterimol$`Only substituent?`,
                                           input_file$Sterimol$`Drop any atoms?`))
    results <- c(results, sterimol.result)
  }

  ### Charges

  if ('Charges' %in% names(input_file)) {
    
    if (length(list.files(pattern = 'nbo.csv',
                          list.dirs(recursive = F))) > 0 &
                          'NPA' %in% input_file$Charges$`Charge Methods`) {
      nbo.result <- list(nbo.df(input_file$Charges$`Atom indices`,
                                paste(input_file$Charges$Differences, collapse = ' ')))
      results <- c(results, nbo.result)
    }

    if (length(list.files(pattern = 'Hirshfeld.csv',
                          list.dirs(recursive = F))) > 0 &
        'Hirshfeld' %in% input_file$Charges$`Charge Methods`) {
      hirsh.result <- list(hirsh.df(input_file$Charges$`Atom indices`,
                                    paste(input_file$Charges$Differences, collapse = ' ')))
      results <- c(results, hirsh.result)
    }
    
    if (length(list.files(pattern = 'CM5.csv',
                          list.dirs(recursive = F))) > 0 &
        'Hirshfeld' %in% input_file$Charges$`Charge Methods`) {
      cm5.result <- list(cm5.df(input_file$Charges$`Atom indices`,
                                  paste(input_file$Charges$Differences, collapse = ' ')))
      results <- c(results, cm5.result)
    }
  }

  ### Dipole

  if ("Dipole Moment" %in% names(input_file)) {
    if (length(input_file$`Dipole Moment`) > 0) {
      dipole.result <- list(dip.gaussian.multi(input_file$`Dipole Moment`))
    } else {
      dipole.result <- list(dip.gaussian.multi())
    }
    results <- c(results, dipole.result)
  }

  ###### Vibrations

  ### Bond Vibrations
  if ("Vibrations" %in% names(input_file)) {
    if ("Bond Vibrations" %in% names(input_file$Vibrations)) {
      bond.vibrations.result <- list(stretch.vib.df(input_file$Vibrations$`Bond Vibrations`))
      results <- c(results, bond.vibrations.result)
    }
  
    ### Ring Vibrations
  
    if ("Ring Vibrations" %in% names(input_file$Vibrations) &
        isTRUE(input_file$Vibrations$`Ring Vibrations`)) {
      ring.vibrations.result <- list(ring.vib.multi(input_file$Vibrations$`Ring Vibrations`))
      results <- c(results, ring.vibrations.result)
    }
  
    ### Bending Vibrations
  
    if ("Bend Vibrations" %in% names(input_file$Vibrations) &
        isTRUE(input_file$Vibrations$`Bend Vibrations`)) {
      bend.vibrations.result <- list(bend.vib.df(input_file$Vibrations$`Bend Vibrations`))
      results <- c(results, bend.vibrations.result)
    }
  }

  ### Angles

  if ("Angles" %in% names(input_file)) {
    angles.result <- list(mol.angles.multi(input_file$Angles))
    results <- c(results, angles.result)
  }

  ### Distances

  if ("Distances" %in% names(input_file)) {
    distances.result <- list(atoms.distance.df(input_file$Distances))
    results <- c(results, distances.result)
  }

  # ### Polarizability
  # 
  # if (input_file$Polarizability[1] == 'yes') {
  #   polar.result <- list(polar.df())
  #   results <- c(results, polar.result)
  # } 



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
Please choose where to save the feature's set file (.csv)",
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


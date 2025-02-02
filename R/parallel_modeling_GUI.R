#' Interactive GUI for Feature Selection and Model Reporting
#'
#' Creates a Shiny web application that provides an interactive interface for feature
#' selection and model reporting. The GUI allows users to upload data, select parameters
#' for model analysis, and generate detailed reports.
#'
#' @title Parallel Model Report GUI
#' @description Launch an interactive Shiny application for model selection and reporting
#'
#' @return A Shiny application object
#' @export
#'
#' @details
#' The GUI provides the following functionality:
#' \itemize{
#'   \item Data upload and preprocessing
#'   \item Sample selection for training/testing
#'   \item Feature selection parameters configuration
#'   \item Cross-validation settings
#'   \item Interactive model selection
#'   \item Results visualization and export
#' }
#'
#' @note 
#' This function requires an active R session with network capabilities for the Shiny
#' interface to work.
#'
#' @importFrom data.table fread rbindlist
#' @importFrom plyr mapvalues
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot geom_point theme_minimal
#' @importFrom ggrepel geom_text_repel
#' @importFrom DT renderDataTable datatable
#' @importFrom tools file_path_sans_ext
model.report.parallel.GUI <- function() {
  # Define UI
  ui <- shiny::fluidPage(
    shiny::titlePanel("Feature Selection and Model Report - Parallel Version"),
    
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::fileInput("data_file", "Upload Data File (CSV)", accept = ".csv"),
        shiny::uiOutput("out_col_ui"),
        
        # Add leave-out method selection
        shiny::radioButtons("leave_out_method", "Leave Out Selection Method:",
                            choices = c("Manual Selection" = "manual",
                                        "Threshold-based" = "threshold")),
        
        # Conditional panels for leave-out methods
        shiny::conditionalPanel(
          condition = "input.leave_out_method == 'manual'",
          shiny::uiOutput("leave_out_ui")
        ),
        shiny::conditionalPanel(
          condition = "input.leave_out_method == 'threshold'",
          shiny::numericInput("leave_out_threshold", "Leave Out Threshold Value", 
                              value = 0, step = 0.1),
          shiny::radioButtons("threshold_direction", "Threshold Direction:",
                              choices = c("Below Threshold" = "below",
                                          "Above Threshold" = "above")),
          shiny::verbatimTextOutput("threshold_preview")
        ),
        
        shiny::numericInput("min", "Minimum # of Variables", value = 2, min = 1),
        shiny::numericInput("max", "Maximum # of Variables", value = NULL),
        shiny::numericInput("cor_threshold", "Correlation Threshold", 
                            value = 0.7, min = 0, max = 1, step = 0.05),
        shiny::uiOutput("folds_ui"),
        shiny::numericInput("iterations", "Number of CV Iterations", value = 1, min = 1),
        shiny::numericInput("cutoff", "Cutoff Value for R^2", value = 0.85, step = 0.05),
        shiny::actionButton("run_model", "Run Model Subset"),
        shiny::actionButton("stop_model", "Stop Execution"),
        shiny::downloadButton("download_results", "Download Model Results (.csv)"),
        shiny::downloadButton("download_included_data", "Download Included Samples (.csv)"),
        shiny::downloadButton("download_excluded_data", "Download Excluded Samples (.csv)"),
        
        # Progress display
        shiny::conditionalPanel(
          condition = "input.run_model > 0",
          shiny::tags$div(
            class = "well",
            shiny::textOutput("progress_status"),
            shiny::div(
              class = "progress",
              shiny::div(id = "modelProgress",
                         class = "progress-bar",
                         role = "progressbar",
                         style = "width: 0%",
                         "aria-valuenow" = "0",
                         "aria-valuemin" = "0",
                         "aria-valuemax" = "100",
                         "0%"
              )
            )
          )
        )
      ),
      
      shiny::mainPanel(
        shiny::tabsetPanel(
          id = "tabs",
          shiny::tabPanel("Model Subset Selection",
                          shiny::uiOutput("metadata_display"),
                          DT::dataTableOutput("results_table")),
          shiny::tabPanel("Model Report",
                          shiny::fileInput("report_model_list", "Upload Results File (CSV)", 
                                           accept = ".csv"),
                          shiny::uiOutput("model_number"),
                          shiny::checkboxInput("predict", "Predict Left-Out Samples", 
                                               value = FALSE),
                          shiny::actionButton("generate_report", "Generate Report"),
                          shiny::plotOutput("model_plot"),
                          shiny::verbatimTextOutput("report_console"))
        )
      )
    )
  )
  
  # Define server
  server <- function(input, output, session) {
    # Progress tracking
    progress <- shiny::reactiveVal(0)
    status_message <- shiny::reactiveVal("")
    
    # Reactive value for stopping execution
    stop_execution <- shiny::reactiveVal(FALSE)
    
    # Reactive value for full data
    full_data_reactive <- shiny::reactive({
      shiny::req(input$data_file)
      data.frame(data.table::fread(input$data_file$datapath), check.names = FALSE)
    })
    
    # Output column selection UI
    output$out_col_ui <- shiny::renderUI({
      shiny::req(full_data_reactive())
      choices <- names(full_data_reactive())[-1]  # Exclude first column
      shiny::selectInput("out_col", "Select Output Column", choices = choices)
    })
    
    # Leave out samples UI
    output$leave_out_ui <- shiny::renderUI({
      shiny::req(full_data_reactive())
      sample_names <- as.character(full_data_reactive()[,1])
      shiny::selectInput("leave_out_manual", "Leave Out Samples", 
                         choices = sample_names, 
                         multiple = TRUE)
    })
    
    # Preview of threshold-based selection
    output$threshold_preview <- shiny::renderPrint({
      shiny::req(selected_leave_out())
      cat("Selected samples to leave out:\n")
      cat(paste(selected_leave_out(), collapse = ", "))
      cat("\n\nTotal samples selected:", length(selected_leave_out()))
    })
    
    # Reactive value for selected leave-out samples
    selected_leave_out <- shiny::reactive({
      shiny::req(full_data_reactive(), input$leave_out_method)
      
      if (input$leave_out_method == "manual") {
        return(input$leave_out_manual)
      } else {
        shiny::req(input$leave_out_threshold, input$threshold_direction, input$out_col)
        df <- full_data_reactive()
        out_values <- as.numeric(df[[input$out_col]])
        sample_names <- df[[1]]
        
        if (input$threshold_direction == "below") {
          return(sample_names[out_values < input$leave_out_threshold])
        } else {
          return(sample_names[out_values > input$leave_out_threshold])
        }
      }
    })
    
    # Reactive value for processed data
    data_reactive <- shiny::reactive({
      shiny::req(full_data_reactive(), input$out_col)
      
      mod_data <- full_data_reactive()
      
      # Save sample names and remove first column
      RN <- as.character(mod_data[,1])
      mod_data <- mod_data[,-1]
      
      # Remove leave-out samples if any
      if (length(selected_leave_out()) > 0) {
        mod_data <- mod_data[!RN %in% selected_leave_out(), ]
        RN <- RN[!RN %in% selected_leave_out()]
      }
      
      # Get column names
      CN <- names(mod_data)
      out.col <- which(CN == input$out_col)
      
      # Scale the data except for the output column
      scaled_data <- data.frame(
        cbind(scale(mod_data[,-out.col], TRUE, TRUE), mod_data[, out.col])
      )
      names(scaled_data)[1:(ncol(scaled_data) - 1)] <- CN[-out.col]
      names(scaled_data)[ncol(scaled_data)] <- CN[out.col]
      row.names(scaled_data) <- RN
      
      # Update max variables input
      shiny::updateNumericInput(session, "max", value = floor(nrow(scaled_data) / 5))
      
      scaled_data
    })
    
    # Reactive value for excluded data
    excluded_data_reactive <- shiny::reactive({
      shiny::req(full_data_reactive(), selected_leave_out())
      
      if (length(selected_leave_out()) > 0) {
        full_data <- full_data_reactive()
        full_data[full_data[,1] %in% selected_leave_out(), ]
      } else {
        NULL
      }
    })
    
    # Number of folds UI
    output$folds_ui <- shiny::renderUI({
      shiny::req(data_reactive())
      remaining_samples <- nrow(data_reactive())
      choices <- c("3-fold" = 3, "5-fold" = 5, "10-fold" = 10)
      
      # Only add LOO if there are not too many samples
      if (remaining_samples <= 100) {
        choices <- c(choices, "LOO" = remaining_samples)
      }
      
      shiny::selectInput("folds", "Number of CV Folds", choices = choices)
    })
    
    # Stop execution handler
    shiny::observeEvent(input$stop_model, {
      stop_execution(TRUE)
    })
    
    # Reset stop execution on new run
    shiny::observeEvent(input$run_model, {
      stop_execution(FALSE)
    })
    
    # Results and metadata reactive
    results_reactive <- shiny::eventReactive(input$run_model, {
      shiny::req(data_reactive())
      
      # Reset progress
      progress(0)
      status_message("Initializing...")
      
      # Create metadata list
      metadata <- list(
        timestamp = Sys.time(),
        data_file = input$data_file$name,
        output_column = input$out_col,
        min_vars = input$min,
        max_vars = input$max,
        folds = as.numeric(input$folds),
        iterations = input$iterations,
        cutoff = input$cutoff,
        correlation_threshold = input$cor_threshold,
        leave_out_method = input$leave_out_method,
        leave_out_threshold = if(input$leave_out_method == "threshold") input$leave_out_threshold else NA,
        threshold_direction = if(input$leave_out_method == "threshold") input$threshold_direction else NA,
        left_out_samples = if(length(selected_leave_out()) > 0) selected_leave_out() else "none"
      )
      
      # Run model subset with error handling
      results <- tryCatch({
        status_message("Running model selection...")
        progress(0.3)
        
        result <- model.subset.parallel(
          data = data_reactive(),
          out.col = which(names(data_reactive()) == input$out_col),
          min = input$min,
          max = input$max,
          results_name = paste0('results_df_', format(Sys.time(), "%Y%m%d_%H%M%S")),
          folds = as.numeric(input$folds),
          iterations = input$iterations,
          cutoff = input$cutoff,
          cor.threshold = input$cor_threshold
        )
        
        result[, 2:4] <- round(result[, 2:4], 2)
        result
        
      }, error = function(e) {
        shiny::showNotification(paste("Error in model processing:", e$message), type = "error")
        NULL
      })
      
      list(
        results = results,
        metadata = metadata
      )
    })
    
    # Results table
    output$results_table <- DT::renderDataTable({
      shiny::req(results_reactive())
      if (stop_execution()) return(NULL)
      
      results <- results_reactive()$results
      if (is.null(results)) return(data.frame(Message = "No results available"))
      
      DT::datatable(
        results,
        options = list(
          pageLength = 15,
          scrollX = TRUE
        ),
        rownames = FALSE
      )
    })
    
    # Metadata display
    output$metadata_display <- shiny::renderUI({
      shiny::req(results_reactive())
      metadata <- results_reactive()$metadata
      
      shiny::div(
        class = "well",
        shiny::h4("Analysis Parameters:"),
        shiny::tags$ul(
          shiny::tags$li(shiny::strong("Timestamp: "), format(metadata$timestamp, "%Y-%m-%d %H:%M:%S")),
          shiny::tags$li(shiny::strong("Data File: "), metadata$data_file),
          shiny::tags$li(shiny::strong("Output Column: "), metadata$output_column),
          shiny::tags$li(shiny::strong("Variables Range: "), paste(metadata$min_vars, "to", metadata$max_vars)),
          shiny::tags$li(shiny::strong("Cross-validation: "), paste(metadata$folds, "folds,", metadata$iterations, "iterations")),
          shiny::tags$li(shiny::strong("R² Cutoff: "), metadata$cutoff),
          shiny::tags$li(shiny::strong("Correlation Threshold: "), metadata$correlation_threshold),
          shiny::tags$li(shiny::strong("Left-out Method: "), metadata$leave_out_method),
          shiny::tags$li(shiny::strong("Left-out Samples: "), paste(metadata$left_out_samples, collapse = ", "))
        )
      )
    })
    
    # Progress outputs
    output$progress_status <- shiny::renderText({
      status_message()
    })
    
    output$modelProgress <- shiny::renderUI({
      shiny::div(
        class = "progress",
        shiny::div(
          class = "progress-bar",
          style = sprintf("width: %d%%;", progress() * 100),
          sprintf("%d%%", progress() * 100)
        )
      )
    })
    
    # Download handlers
    output$download_results <- shiny::downloadHandler(
      filename = function() {
        paste0("model_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        utils::write.csv(results_reactive()$results, file, row.names = FALSE)
      },
      contentType = "text/csv"
    )
    
    output$download_included_data <- shiny::downloadHandler(
      filename = function() {
        paste0("included_samples_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        utils::write.csv(data_reactive(), file)
      },
      contentType = "text/csv"
    )
    
    output$download_excluded_data <- shiny::downloadHandler(
      filename = function() {
        paste0("excluded_samples_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        excluded_data <- excluded_data_reactive()
        if (!is.null(excluded_data)) {
          utils::write.csv(excluded_data, file)
        } else {
          utils::write.csv(data.frame(Message = "No samples were excluded"), file)
        }
      },
      contentType = "text/csv"
    )
    
    # Model report generation
    output$model_number <- shiny::renderUI({
      shiny::req(input$report_model_list)
      model_data <- data.frame(data.table::fread(input$report_model_list$datapath))
      shiny::selectInput("model.num", "Select Model Number",
                         choices = seq_len(nrow(model_data)),
                         selected = 1)
    })
    
    model_report_reactive <- shiny::eventReactive(input$generate_report, {
      shiny::req(input$data_file, input$report_model_list, input$model.num)
      
      tryCatch({
        model.report.from.list(
          dataset = input$data_file$datapath,
          model.list = input$report_model_list$datapath,
          out.col = input$out_col,
          leave.out = if(length(selected_leave_out()) > 0) selected_leave_out() else '',
          predict = input$predict,
          what.model = as.numeric(input$model.num)
        )
      }, error = function(e) {
        shiny::showNotification(paste("Error in model report generation:", e$message), type = "error")
        NULL
      })
    })
    
    output$model_plot <- shiny::renderPlot({
      shiny::req(model_report_reactive())
      model_report_reactive()
    })
    
    output$report_console <- shiny::renderPrint({
      shiny::req(model_report_reactive())
    })
  }
  
  # Create and run the Shiny app
  shiny::shinyApp(ui = ui, server = server)
}

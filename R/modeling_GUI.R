##### Linear modeling GUI - allows users to work with data sets and retrieve 
##### linear regression reports

#' Screen for models, cross validate them and plot. Designed for interactive work.
#' @return csv file with model list, and a plot
#' @importFrom utils capture.output
#' @export
model.report.GUI <- function() {
  # Define UI for the app
  ui <- shiny::fluidPage(
    shiny::titlePanel("Feature Selection and Model Report"),
    
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::fileInput("data_file", "Upload Data File (CSV)", accept = ".csv"),
        shiny::uiOutput("out_col_ui"),  # Dynamic output column selection
        shiny::uiOutput("leave_out_ui"),  # Dynamic leave out samples dropdown
        shiny::numericInput("min", "Minimum # of Variables", value = 2, min = 1),
        shiny::numericInput("max", "Maximum # of Variables", value = NULL),
        shiny::uiOutput("folds_ui"),
        shiny::numericInput("iterations", "Number of CV Iterations", value = 1, min = 1),
        shiny::numericInput("cutoff", "Cutoff Value for R^2", value = 0.85, step = 0.05),
        shiny::actionButton("run_model", "Run Model Subset"),
        shiny::actionButton("stop_model", "Stop Execution"),
        shiny::downloadButton("download_results", "Download Results as CSV")
      ),
      
      shiny::mainPanel(
        shiny::tabsetPanel(
          id = "tabs",
          shiny::tabPanel("Model Subset Selection",
                          shiny::tableOutput("results_table")),
          shiny::tabPanel("Model Report",
                          shiny::fileInput("report_model_list", "Upload Results File (CSV)", accept = ".csv"),
                          shiny::uiOutput("model_number"),
                          shiny::checkboxInput("predict", "Predict Left-Out Samples", value = FALSE),
                          shiny::actionButton("generate_report", "Generate Report"),
                          shiny::plotOutput("model_plot"),
                          shiny::verbatimTextOutput("report_console"))
        )
      )
    )
  )
  
  # Define server logic
  server <- function(input, output, session) {
    # Reactive value to track if the process should be stopped
    stop_execution <- shiny::reactiveVal(FALSE)
    
    data_reactive <- shiny::reactive({
      shiny::req(input$data_file)
      mod_data <- data.frame(data.table::fread(input$data_file$datapath), check.names = F)
      
      # Preprocessing steps
      RN <- mod_data[, 1]  # Save row names
      mod_data <- mod_data[, -1]  # Remove the first column
      mod_data <- mod_data[complete.cases(mod_data), ]  # Remove rows with NA values
      CN <- names(mod_data)  # Save column names
      shiny::req(input$out_col)  # Ensure output column is selected
      out.col <- which(CN == input$out_col)
      # Scale the data except for the output column
      mod_data <- data.frame(cbind(scale(mod_data[,-out.col], T, T), mod_data[, out.col]))
      names(mod_data)[1:(ncol(mod_data) - 1)] <- CN[-out.col]
      names(mod_data)[ncol(mod_data)] <- CN[out.col]  # Restore column names
      row.names(mod_data) <- RN  # Restore row names
      
      # Update dynamic default values after loading data
      shiny::updateNumericInput(session, "max", value = floor(nrow(mod_data) / 5))
      
      mod_data
    })
    
    # Dynamic output column selection
    output$out_col_ui <- shiny::renderUI({
      shiny::req(input$data_file)
      shiny::selectInput("out_col", "Select Output Column", 
                         choices = names(data.frame(data.table::fread(input$data_file$datapath), check.names = F))[-1])
    })
    
    # Dynamic leave out samples selection
    output$leave_out_ui <- shiny::renderUI({
      shiny::req(data_reactive())
      shiny::selectInput("leave_out", "Leave Out Samples", 
                         choices = data.frame(data.table::fread(input$data_file$datapath), check.names = F)[, 1], 
                         multiple = TRUE)
    })
    
    # Dynamic number of folds selection
    output$folds_ui <- shiny::renderUI({
      shiny::req(data_reactive())
      shiny::selectInput("folds", "Number of CV Folds", 
                         choices = c("3-fold" = 3, "5-fold" = 5, "10-fold" = 10, 
                                     "LOO" = nrow(data_reactive())))
    })
    
    shiny::observeEvent(input$run_model, {
      stop_execution(FALSE)  # Reset the stop flag when the model is run
    })
    
    results_reactive <- shiny::eventReactive(input$run_model, {
      shiny::req(data_reactive())
      
      model.subset(data = data_reactive(),
                   out.col = which(names(data_reactive()) == input$out_col),
                   min = input$min,
                   max = input$max,
                   folds = as.numeric(input$folds),
                   iterations = input$iterations,
                   cutoff = input$cutoff,
                   leave.out = input$leave_out)
    })
    
    shiny::observeEvent(input$stop_model, {
      stop_execution(TRUE)
    })
    
    output$results_table <- shiny::renderTable({
      shiny::req(results_reactive())
      if (stop_execution()) {
        return(NULL)
      }
      results_reactive()
    })
    
    output$download_results <- shiny::downloadHandler(
      filename = function() {
        paste("model_subset_results", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        utils::write.csv(results_reactive(), file, row.names = FALSE)
      }
    )
    
    # Model report section - Fixed version
    model_list_reactive <- shiny::reactive({
      shiny::req(input$report_model_list)
      model_list <- data.table::fread(input$report_model_list$datapath)
      # Ensure the Model column exists and contains valid indices
      if (!"Model" %in% names(model_list)) {
        model_list$Model <- seq_len(nrow(model_list))
      }
      model_list
    })
    
    output$model_number <- shiny::renderUI({
      shiny::req(model_list_reactive())
      models <- model_list_reactive()$Model
      shiny::selectInput("model.num", "Select Model Number", 
                         choices = models,
                         selected = models[1])
    })
    
    # Create a reactive value to store the report results
    report_results <- shiny::reactiveVal(NULL)
    
    shiny::observeEvent(input$generate_report, {
      shiny::req(input$data_file, input$report_model_list, input$model.num)
      
      # Ensure the data is properly loaded
      data <- data.frame(data.table::fread(input$data_file$datapath), check.names = F)
      model_list <- data.table::fread(input$report_model_list$datapath)
      
      # Generate the report
      tryCatch({
        # Debug prints
        print("Debug info:")
        print(paste("Dataset path:", input$data_file$datapath))
        print(paste("Model list path:", input$report_model_list$datapath))
        print(paste("Output column:", input$out_col))
        print(paste("Leave out:", paste(input$leave_out, collapse=", ")))
        print(paste("Model number:", input$model.num))
        
        result <- model.report.from.list(
          dataset = input$data_file$datapath,  # Pass the file path instead of data frame
          model.list = input$report_model_list$datapath,  # Pass the file path instead of data frame
          out.col = input$out_col,
          leave.out = input$leave_out,
          predict = input$predict,
          what.model = as.numeric(input$model.num)
        )
        
        report_results(result)
        
      }, error = function(e) {
        # Handle any errors that occur during report generation
        shiny::showNotification(
          paste("Error generating report:", e$message),
          type = "error"
        )
      })
    })
    
    # Render the plot
    output$model_plot <- shiny::renderPlot({
      shiny::req(report_results())
      if (inherits(report_results(), "ggplot")) {
        report_results()
      } else if (is.list(report_results()) && "plot" %in% names(report_results())) {
        report_results()$plot
      }
    })
    
    # Render the console output
    output$report_console <- shiny::renderPrint({
      shiny::req(report_results())
      if (is.list(report_results()) && "summary" %in% names(report_results())) {
        report_results()$summary
      } else {
        utils::capture.output(report_results())
      }
    })
  }
  
  # Run the application 
  shiny::shinyApp(ui = ui, server = server)
}
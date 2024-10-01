##### Linear modeling GUI - allows users to work with data sets and retrieve 
##### linear regression reports

#' Screen for models, cross validate them and plot. Designed for interactive work.
#' @return csv file with model list, and a plot
#'  use and documentation
#' @export
model.report.GUI <- function() {
  # Define UI for the app
  ui <- fluidPage(
    titlePanel("Feature Selection and Model Report"),
    
    sidebarLayout(
      sidebarPanel(
        fileInput("data_file", "Upload Data File (CSV)", accept = ".csv"),
        uiOutput("out_col_ui"),  # Dynamic output column selection
        uiOutput("leave_out_ui"),  # Dynamic leave out samples dropdown
        numericInput("min", "Minimum # of Variables", value = 2, min = 1),
        numericInput("max", "Maximum # of Variables", value = NULL),  # Restored maximum variables input
        uiOutput("folds_ui"),  # Dynamic number of folds dropdown
        numericInput("iterations", "Number of CV Iterations", value = 1, min = 1),
        numericInput("cutoff", "Cutoff Value for R^2", value = 0.85, step = 0.05),
        actionButton("run_model", "Run Model Subset"),
        actionButton("stop_model", "Stop Execution"),  # Added stop execution button
        downloadButton("download_results", "Download Results as CSV")  # Added download button
      ),
      
      mainPanel(
        tabsetPanel(
          id = "tabs",
          tabPanel("Model Subset Selection",
                   tableOutput("results_table")),
          tabPanel("Model Report",
                   fileInput("report_model_list", "Upload Results File (CSV)", accept = ".csv"),
                   uiOutput("model_number"),
                   checkboxInput("predict", "Predict Left-Out Samples", value = FALSE),
                   actionButton("generate_report", "Generate Report"),
                   plotOutput("model_plot"),
                   verbatimTextOutput("report_console"))
        )
      )
    )
  )
  
  # Define server logic
  server <- function(input, output, session) {
    # Reactive value to track if the process should be stopped
    stop_execution <- reactiveVal(FALSE)
    
    data_reactive <- reactive({
      req(input$data_file)
      mod_data <- data.frame(data.table::fread((input$data_file$datapath)),check.names = F)
      
      # Preprocessing steps
      RN <- mod_data[, 1]  # Save row names
      mod_data <- mod_data[, -1]  # Remove the first column
      mod_data <- mod_data[complete.cases(mod_data), ]  # Remove rows with NA values
      CN <- names(mod_data)  # Save column names
      out.col <- which(CN == input$out_col)
      # Scale the data except for the output column
      mod_data <- data.frame(cbind(scale(mod_data[,-out.col], T, T), mod_data[, out.col]))
      names(mod_data)[1:(ncol(mod_data) - 1)] <- CN[-out.col]
      names(mod_data)[ncol(mod_data)] <- CN[out.col]  # Restore column names
      row.names(mod_data) <- RN  # Restore row names
      
      # Update dynamic default values after loading data
      updateNumericInput(session, "max", value = floor(nrow(mod_data) / 5))
      
      mod_data
    })
    
    # Dynamic output column selection
    output$out_col_ui <- renderUI({
      req(input$data_file)
      selectInput("out_col", "Select Output Column", 
                  choices = names(data.frame(data.table::fread((input$data_file$datapath)),check.names = F))[-1])
    })
    
    # Dynamic leave out samples selection
    output$leave_out_ui <- renderUI({
      req(data_reactive())
      selectInput("leave_out", "Leave Out Samples", 
                  choices = data.frame(data.table::fread((input$data_file$datapath)), check.names = F)[ ,1], multiple = TRUE)
    })
    
    # Dynamic number of folds selection
    output$folds_ui <- renderUI({
      selectInput("folds", "Number of CV Folds", 
                  choices = c("3-fold" = 3, "5-fold" = 5, "10-fold" = 10, "LOO" = nrow(data_reactive())))
    })
    
    observeEvent(input$run_model, {
      stop_execution(FALSE)  # Reset the stop flag when the model is run
    })
    
    results_reactive <- eventReactive(input$run_model, {
      req(data_reactive())
      
      model.subset(data = data_reactive(),
                   out.col = which(names(data_reactive()) == input$out_col),  # Get the index of the selected column
                   min = input$min,
                   max = input$max,
                   folds = as.numeric(input$folds),  # Convert folds selection to numeric
                   iterations = input$iterations,
                   cutoff = input$cutoff,
                   leave.out = input$leave_out)  # Add leave.out argument
    })
    
    observeEvent(input$stop_model, {
      stop_execution(TRUE)  # Set the stop flag when the stop button is pressed
    })
    
    output$results_table <- renderTable({
      req(results_reactive())
      if (stop_execution()) {
        return(NULL)  # Stop execution if the stop flag is set
      }
      results_reactive()
    })
    
    # Download handler for saving results as a CSV file
    output$download_results <- downloadHandler(
      filename = function() {
        paste("model_subset_results", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(results_reactive(), file, row.names = FALSE)
      }
    )
    
    # Model report section
    
    
    model_list_reactive <- reactive({
      req(input$report_model_list)
      data.table::fread(input$report_model_list$datapath)
    })
    
    output$model_number <- renderUI({
      req(input$report_model_list)
      selectInput("model.num", "Select what Model to Report", choices = row.names(model_list_reactive()), multiple = FALSE)
    })
    
    observeEvent(input$generate_report, {
      req(model_list_reactive())
      
      # Validate inputs
      req(input$data_file$datapath, input$report_model_list)
      
      # Call the model report function
      result <- model.report.from.list(
        dataset = input$data_file$datapath,
        model.list = input$report_model_list$datapath,
        out.col = input$out_col,
        leave.out = input$leave_out,  # Use leave.out from the first section
        predict = input$predict,
        what.model = input$model.num
      )
      
      # Update plot and console
      output$model_plot <- renderPlot({
        result
      })
      
      output$report_console <- renderPrint({
        capture.output(result)
      })
    })
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
}

####### ----------------------------------------------------#####
####### -----------------Utility Functions------------------#####
####### ----------------------------------------------------#####

#' Cross validate (k-fold) a single model
#'
#' A utility function for the brute force search of models
#' @param formula model formula (not designed for direct use)
#' @param data data frame with outcome column
#' @param out.col column number with output
#' @param folds define K
#' @keywords internal
#' @return MAE and Q2
model.single.cv <- function(formula, data, out.col, folds) {
  predictions <- list()
  if (folds == nrow(data)) {
    split.assign <- sample(1:folds, nrow(data), replace = F)
  } else {
    split.assign <- sample(1:folds, nrow(data), replace = T)
  }
  cvdat <- cbind(data, split.assign)
  colnames(cvdat)[dim(data)[2] + 1] <- "assign"
  for (i in 1:folds) {
    train <- cvdat[cvdat$assign != i, ]
    test <- cvdat[cvdat$assign == i, ]
    model <- lm(formula, data = train)
    predictions[[i]] <- data.frame(predict(model, newdata = test))
  }
  pred <- dplyr::bind_rows(predictions)
  pred <- pred[match(rownames(cvdat), rownames(pred)), ]
  sq.dis.cv <- data.frame(abs(pred - data[, out.col]))
  MAE.3cv <- mean(sq.dis.cv[, 1])
  q2.3cv <- caret::R2(pred, data[, out.col])
  return(list(MAE.3cv, q2.3cv))
}

#' Cross validate (k-fold) a single model
#'
#' An iterative version of model.single.cv
#' @param formula model formula (not designed for direct use)
#' @param data data frame with outcome column
#' @param out.col column number with output
#' @param folds define K
#' @param iterations Number of iterations
#' @export
#' @return averaged MAE and Q2
model.cv <- function(formula, data, out.col, folds, iterations) {
  mae.list <- list()
  q2.list <- list()
  for (i in 1:iterations) {
    tool <- model.single.cv(formula, data, out.col, folds)
    mae.list[[i]] <- tool[[1]]
    q2.list[[i]] <- tool[[2]]
  }
  MAE.validation <- Reduce(`+`, mae.list) / iterations
  q2.validation <- Reduce(`+`, q2.list) / iterations
  return(list(MAE.validation, q2.validation))
}

#' Create Subsets of Features Based on Correlation Threshold and R-squared
#'
#' This function generates all possible subsets of features in a user-defined range (min to max), while ensuring
#' no two variables in a subset have a Pearson correlation coefficient of `cor.threshold` or higher. 
#' It returns the top models based on R-squared and further evaluates them using cross-validation to calculate Q-squared and MAE.
#'
#' @param data A data frame containing the dataset.
#' @param leave.out A vector of row names to exclude from the data (default is an empty string).
#' @param out.col The index of the output column in the data (default is the last column).
#' @param min Minimum number of variables in a subset (default is 2).
#' @param max Maximum number of variables in a subset (default is floor(n/5)).
#' @param folds Number of folds for cross-validation (default is number of rows in data).
#' @param iterations Number of iterations for cross-validation (default is 1).
#' @param cutoff R-squared cutoff for model selection (default is 0.85).
#' @param cor.threshold Pearson correlation threshold to filter highly correlated variables (default is 0.7).
#' @export
#' @return A data frame containing the top 10 models based on R-squared, Q-squared, and MAE.
model.subset <- function(data, leave.out = '', out.col = dim(data)[2],
                         min = 2, max = floor(dim(data)[1] / 5),
                         folds = nrow(data), iterations = 1,
                         cutoff = 0.85, cor.threshold = 1) {
  
  # Exclude specified rows
  data <- data[!(row.names(data) %in% leave.out), ]
  
  # Get output variable name
  output <- stringr::str_c("`", names(data[out.col]), "`")
  
  # Get names of all predictor variables
  vars <- names(data[, -out.col])
  for (i in 1:length(vars)) {
    vars[i] <- stringr::str_c("`", vars[i], "`")
  }
  
  comb.list <- list()
  ols.list <- list()
  q2.list <- list()
  mae.list <- list()
  
  # Calculate correlation matrix and apply the correlation threshold filter
  cor.matrix <- cor(data[, -out.col], method = "pearson")
  
  for (i in min:max) {
    # Get combinations of variables for the current subset size
    combinations <- aperm(combn(vars, i))
    valid.combinations <- list()
    
    # Check correlation threshold for each combination
    for (j in 1:nrow(combinations)) {
      subset.vars <- gsub("`", "", combinations[j, ])
      
      # Extract correlation values for the subset of variables
      cor.subset <- cor.matrix[subset.vars, subset.vars]
      cor.upper <- cor.subset[upper.tri(cor.subset)]
      
      # If all correlations are below the threshold, add to valid combinations
      if (all(abs(cor.upper) < cor.threshold)) {
        valid.combinations[[length(valid.combinations) + 1]] <- combinations[j, ]
      }
    }
    
    # Store valid combinations
    comb.list <- do.call(rbind, valid.combinations)
    
    # Create formulas for valid combinations
    for (i in 1:max) {
      comb.list[[i]] <- data.frame(aperm(combn(vars, i)), stringsAsFactors = F)
      comb.list[[i]][, dim(comb.list[[i]])[2] + 1] <- do.call(
        paste,
        c(comb.list[[i]][names(comb.list[[i]])],
          sep = " + "
        )
      )
      names(comb.list[[i]])[dim(comb.list[[i]])[2]] <- "formula"
      for (co in names(comb.list[[i]])[1:length(names(comb.list[[i]])) - 1]) {
        comb.list[[i]][co] <- NULL
      }
    }
    comb.list <- plyr::compact(comb.list)
    forms <- do.call(rbind, comb.list)
    forms$formula <- stringr::str_c(output, " ~ ", forms$formula)
    
    # Evaluate R-squared for each formula
    ols.list <- lapply(forms$formula, function(x) summary(lm(x, data = data))$r.squared)
    forms[, 2] <- do.call(rbind, ols.list)
    names(forms)[2] <- "R.sq"
    
    # Apply R-squared cutoff
    forms.cut <- forms[forms$R.sq > cutoff, ]
    
    # If not enough models pass the cutoff, lower the cutoff incrementally
    while (nrow(forms.cut) <= 10) {
      cutoff <- cutoff - 0.1
      forms.cut <- forms[forms$R.sq > cutoff, ]
    }
    
    # Select the top 10 models by R-squared
    forms.cut <- dplyr::arrange(forms.cut, desc(forms.cut$R.sq))
    if (nrow(forms.cut) >= 10) forms.cut <- forms.cut[1:10, ]
    
    # Perform cross-validation for each model to calculate Q-squared and MAE
    for (i in 1:nrow(forms.cut)) {
      stts <- model.cv(forms.cut[i, 1], data, out.col, folds, iterations)
      q2.list[[i]] <- stts[2]
      mae.list[[i]] <- stts[1]
    }
    
    # Add Q-squared and MAE to the results
    forms.cut[, 3] <- data.table::transpose(do.call(rbind, q2.list))
    forms.cut[, 4] <- data.table::transpose(do.call(rbind, mae.list))
    names(forms.cut)[3:4] <- c("Q.sq", "MAE")
    
    # Add model numbering and arrange by Q-squared
    forms.cut <- dplyr::arrange(forms.cut, desc(forms.cut$Q.sq))
    forms.cut <- dplyr::mutate(forms.cut, Model = seq(1, nrow(forms.cut)))
  }
  
  return(forms.cut)
}

#' ggplot2 plot of a linear model for QSAR
#'
#' @param model model formula (not designed for direct use)
#' @param data data frame with outcome column
#'
#' @return a chemistry relevant plot
#' @export
model.plot <- function(model = models[1,1], data) {
  best.mod <- lm(model, data = data)
  pred_interval <- predict(best.mod,
                           newdata = mod_data,
                           interval = 'pre',
                           level = 0.9)
  plot.dat <- data.frame(cbind(mod_data[dim(mod_data)[2]], pred_interval))
  colnames(plot.dat) <- c('Measured', 'Predicted', 'lwr', 'upr')
  rownames(plot.dat) <- row.names(mod_data)
  
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o_",'2-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"m_",'3-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"p_",'4-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o4-",'2,4-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"m3-",'3,5-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o3-",'2,3-')
  
  plot.dat <- dplyr::mutate(plot.dat, Position = rep(NA, nrow(plot.dat)))
  
  for (i in 1:nrow(mod_data)) {
    if (grepl('3-',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'meta'
    }
    if (grepl('5-',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'meta'
    }
    if (grepl('2-',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'ortho'
    }
    if (grepl('basic',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'Ph'
    }
    if (grepl('penta_F',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'C6F5'
    }
    if (grepl('4-',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'para'
    }
  }
  plot.dat <- dplyr::mutate(plot.dat, label = row.names(plot.dat))
  
  plot.dat <- dplyr::mutate(plot.dat, 
                            shapes = c(rep(18, nrow(mod_data))))
  
  
  # Calculate absolute errors
  plot.dat$error <- abs(plot.dat$Measured - plot.dat$Predicted)
  # Create label vector that only includes points with large errors
  plot.dat$show_label <- ifelse(plot.dat$error > 15, plot.dat$label, "")
  
  # Extract dataset name (without path and extension) for subtitle
  dataset_name <- basename(dataset)
  dataset_name <- sub("\\.[^.]*$", "", dataset_name)
  
  # Get model formula as string for subtitle
  model_formula_str <- as.character(models$formula[what.model])
  
  # Store for plot subtitles
  plot_subtitle1 <- paste("Model:", model_formula_str)
  plot_subtitle2 <- paste("Dataset:", dataset_name)
  
  # Calculate plot limits for proper legend positioning and for x/y limits
  x_min <- min(plot.dat[1:nrow(mod_data),3])
  x_max <- max(plot.dat[1:nrow(mod_data),4])
  y_min <- min(plot.dat[1:nrow(mod_data),3])
  y_max <- max(plot.dat[1:nrow(mod_data),4])
  
  plot <- suppressMessages(ggplot2::ggplot(plot.dat, ggplot2::aes(x = Measured, y = Predicted)) +
                             ggplot2::geom_point(size = 2, shape = plot.dat$shapes, ggplot2::aes(color = Position)) +
                             ggplot2::stat_smooth(ggplot2::aes(y = lwr), color = "cadetblue", linetype = "dashed",
                                                  se = F, method = 'lm', fullrange = T, size = 0.8) +
                             ggplot2::stat_smooth(ggplot2::aes(y = upr), color = "cadetblue", linetype = "dashed",
                                                  se = F, method = 'lm', fullrange = T, size = 0.8) +
                             ggplot2::labs(x = 'Measured',
                                           y = 'Predicted',
                                           title = plot_title, 
                                           subtitle = paste(plot_subtitle1, "\n", plot_subtitle2)) +
                             ggplot2::stat_smooth(method = 'lm',se = F, formula = y~x,
                                                  color = 'black',fullrange = T, linetype = 'dashed') +
                             ggplot2::theme(axis.line.x = ggplot2::element_line(linewidth = 1, colour = "black"),
                                            axis.line.y = ggplot2::element_line(linewidth = 1, colour = "black"),
                                            axis.text.x = ggplot2::element_text(colour = "black", size = 12,face = 'bold'),
                                            axis.text.y = ggplot2::element_text(colour = "black", size = 12,face = 'bold'),
                                            axis.title.x = ggplot2::element_text(colour = "black", size = 12,face = 'bold'),
                                            axis.title.y = ggplot2::element_text(colour = "black", size = 12,face = 'bold'),
                                            panel.grid.major = ggplot2::element_blank(),
                                            panel.grid.minor = ggplot2::element_blank(),
                                            panel.border = ggplot2::element_blank(), 
                                            panel.background = ggplot2::element_blank(),
                                            # Position legend at the bottom right inside the plot
                                            legend.position = c(0.85, 0.15),
                                            legend.justification = c(0.5, 0.5),
                                            legend.background = ggplot2::element_blank(),
                                            legend.key = ggplot2::element_blank(),
                                            legend.key.size = unit(0.8, "lines"),
                                            legend.title = ggplot2::element_text(size = 10),
                                            legend.text = ggplot2::element_text(size = 9),
                                            # Title and subtitle styling
                                            plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                            plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5)) +
                             ggplot2::scale_color_manual('', values = c(Ph = "black", meta = 'tan1', C6F5 = 'darkgrey',
                                                                        para = '#66a182',ortho = '#d1495b', external = 'steelblue4')) +
                             ggplot2::xlim(x_min, x_max) +
                             ggplot2::ylim(y_min, y_max) +
                             ggplot2::coord_fixed(ratio = 1) +
                             ggrepel::geom_text_repel(
                               data = subset(plot.dat, error > 13),
                               aes(label = label),
                               size = 3,
                               min.segment.length = Inf,
                               seed = 42,
                               point.padding = 0.4,
                               segment.color = 'grey50',
                               force_pull = 0.02,
                               nudge_x = 0.022,
                               direction = 'y'
                             ) +
                             ggplot2::theme(text = ggplot2::element_text(family = 'Arial')) +
                             ggplot2::annotate('text',
                                               x = min(plot.dat[1:nrow(mod_data),3]),
                                               y = max(plot.dat[1:nrow(mod_data),2]), label = annotations,
                                               parse = F,
                                               hjust = "left", vjust = 0))
  plot
}

####### ----------------------------------------------------#####
####### -------------------User Functions-------------------#####
####### ----------------------------------------------------#####

#' Create a Formatted Table of Out-of-Sample Validation Results
#'
#' This function generates a formatted ggplot visualization of out-of-sample prediction results,
#' presenting the predicted values, measured values, and errors in a tabular format. The table 
#' includes color-coded error indicators to highlight which predictions meet an acceptable error threshold.
#'
#' @param oos_data A data frame containing out-of-sample validation results with the following columns:
#'        \itemize{
#'          \item \code{OOS Pred} - Predicted values for out-of-sample data
#'          \item \code{OOS Measured} - Actual measured values for out-of-sample data
#'          \item \code{OOS Error} - Difference between measured and predicted values
#'        }
#' @param plot.title Title for the plot (default: "Out-of-Sample Validation")
#' @param subtitle Subtitle for the plot (default: "Predicted vs Measured Values")
#' @param conformation Optional text to indicate conformational information, added to the caption (default: "")
#' @param error_threshold Numeric threshold for acceptable prediction error (default: 20). Errors with 
#'        absolute values below or equal to this threshold are marked as acceptable.
#'
#' @return A ggplot object representing a formatted table with the following features:
#'         \itemize{
#'           \item Row names from the input data frame as compound identifiers
#'           \item Columns for predicted values, experimental values, and errors
#'           \item Color-coded errors (green for acceptable, red for unacceptable)
#'           \item Summary statistics in the caption (MAE and percentage of predictions within threshold)
#'         }
#'
#' @details
#' The function creates a custom ggplot visualization that mimics a table layout with grid lines.
#' Row heights are equal, and the table is designed with a fixed aspect ratio for consistent display.
#' Errors are color-coded based on the specified threshold to provide a quick visual assessment
#' of prediction quality.
#'
#' The table includes a caption with summary statistics:
#' \itemize{
#'   \item Mean Absolute Error (MAE) across all predictions
#'   \item Percentage of predictions that fall within the acceptable error threshold
#'   \item Optional conformational information if provided
#' }
#'
#' @importFrom ggplot2 ggplot geom_segment geom_text xlim ylim theme_void theme 
#' @importFrom ggplot2 element_text margin labs coord_fixed scale_color_identity aes
#'
#' @export
oos_validation_table <- function(oos_data, 
                                 plot.title = "Out-of-Sample Validation",
                                 subtitle = "Predicted vs Measured Values",
                                 conformation = "",
                                 error_threshold = 20) {  # Threshold to consider errors acceptable
  
  # Required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required")
  }
  
  # Make a copy of the input data
  df <- oos_data
  
  # Add a column to mark errors below threshold as acceptable
  df$Acceptable <- abs(df$`OOS Error`) <= error_threshold
  
  # Calculate the percentage of acceptable predictions
  accuracy_pct <- round(mean(df$Acceptable) * 100, 1)
  
  # Calculate mean absolute error (MAE)
  mae <- round(mean(abs(df$`OOS Error`)), 1)
  
  # Number of rows
  n_rows <- nrow(df)
  
  # Format the ggplot - starting with an empty plot
  plot <- ggplot2::ggplot() +
    
    # Set up the grid for the table - with equal-height cells
    # Bottom of table at y=0, top of table at y=n_rows
    
    # Vertical grid lines (columns)
    ggplot2::geom_segment(
      data = data.frame(x = c(0, 1, 2, 3, 4)),
      mapping = ggplot2::aes(
        x = x, xend = x,
        y = 0, yend = n_rows
      ),
      color = "black", size = 0.5
    ) +
    
    # Horizontal grid lines (rows) - with EQUAL height cells
    ggplot2::geom_segment(
      data = data.frame(y = 0:n_rows),
      mapping = ggplot2::aes(
        x = 0, xend = 4,
        y = y, yend = y
      ),
      color = "black", size = 0.5
    ) +
    
    # Column headers directly at the top of the table
    ggplot2::geom_text(
      data = data.frame(
        x = c(0.5, 1.5, 2.5, 3.5),
        y = rep(n_rows + 0.5, 4),
        label = c("", "Pred", "Exp", "Error"),
        stringsAsFactors = FALSE
      ),
      mapping = ggplot2::aes(
        x = x, y = y,
        label = label
      ),
      size = 4.5, fontface = "bold"
    ) +
    
    # Compound names
    ggplot2::geom_text(
      data = data.frame(
        x = rep(0.5, n_rows),
        y = n_rows:1 - 0.5,  # Reverse order for y to match the table image
        label = rownames(df),
        stringsAsFactors = FALSE
      ),
      mapping = ggplot2::aes(
        x = x, y = y,
        label = label
      ),
      size = 4, fontface = "bold"
    ) +
    
    # Predicted values
    ggplot2::geom_text(
      data = data.frame(
        x = rep(1.5, n_rows),
        y = n_rows:1 - 0.5,  # Reverse order for y
        label = df$`OOS Pred`,
        stringsAsFactors = FALSE
      ),
      mapping = ggplot2::aes(
        x = x, y = y,
        label = label
      ),
      size = 4
    ) +
    
    # Measured values
    ggplot2::geom_text(
      data = data.frame(
        x = rep(2.5, n_rows),
        y = n_rows:1 - 0.5,  # Reverse order for y
        label = df$`OOS Measured`,
        stringsAsFactors = FALSE
      ),
      mapping = ggplot2::aes(
        x = x, y = y,
        label = label
      ),
      size = 4
    ) +
    
    # Error values with color coding
    ggplot2::geom_text(
      data = data.frame(
        x = rep(3.5, n_rows),
        y = n_rows:1 - 0.5,  # Reverse order for y
        label = df$`OOS Error`,
        acceptable = df$Acceptable,
        stringsAsFactors = FALSE
      ),
      mapping = ggplot2::aes(
        x = x, y = y,
        label = label,
        color = ifelse(acceptable, "#008000", "#FF0000")
      ),
      size = 4, fontface = "bold"
    ) +
    
    # Use colors directly
    ggplot2::scale_color_identity() +
    
    # Set limits to include headers above the table
    ggplot2::xlim(-0.5, 4.5) +
    ggplot2::ylim(-0.5, n_rows + 1) +
    
    # Theme customization
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 14, hjust = 0.5),
      # Center the caption with hjust = 0.5 and add top margin
      plot.caption = ggplot2::element_text(size = 10, hjust = 0.5, margin = ggplot2::margin(t = 15)),
      plot.margin = ggplot2::margin(15, 15, 15, 15)
    ) +
    
    # Add title, subtitle and caption
    ggplot2::labs(
      title = plot.title,
      subtitle = subtitle,
      caption = paste0(
        "Mean Absolute Error: ", mae, 
        " | Predictions within threshold: ", accuracy_pct, "%",
        if(nchar(conformation) > 0) paste0("\n", conformation) else ""
      )
    ) +
    
    # Set aspect ratio - make cells rectangular but maintain fixed width-to-height ratio
    ggplot2::coord_fixed(ratio = 0.33, clip = "off")
  
  return(plot)
}

#'  Generate models, cross validate, and create visualizations
#'
#' Screen for models, cross validate them and plot. Designed for interactive work.
#' @param dataset A dataframe with outcome column
#' @param min Minimum # of features (default = 2)
#' @param max Max # of features (defaults = # of observations / 5)
#' @param leave.out Name of observations to leave out (e.g. 'p_Br')
#' @param predict If leave.out is not empty, should a prediction of it be computed?
#' @param what.model If not in an interactive session, what model should be used?
#' @param cor.threshold Pearson correlation threshold for variable selection (default = 0.7)
#' @param save.pred Logical indicating whether to save out-of-sample predictions to a CSV file (default: TRUE)
#' @param plot_title Title for the generated plot (default: "Linear Regression Model Analysis")
#' @return models list, CV results for k=3/5/LOO and a plot of choice. interactive.
#' @export
model.report <- function(dataset, min = 2, max = NULL,
                         leave.out = '', predict = FALSE,
                         what.model = NULL, cor.threshold = 1,
                         save.pred = TRUE,
                         plot_title = "Linear Regression Model Analysis") {
  # Setup default check.names behavior
  default::default(data.frame) <- list(check.names = FALSE)
  
  # Extract dataset name for output information
  dataset_name <- tools::file_path_sans_ext(basename(dataset))
  cat("Dataset:", dataset_name, "\n")
  
  # Read and preprocess data
  tryCatch({
    mod_data <- data.frame(data.table::fread(dataset, header = TRUE, check.names = FALSE))
  }, error = function(e) {
    stop(paste("Error reading dataset:", e$message))
  })
  
  # Extract row names (first column) and remove it from the dataframe
  RN <- mod_data[,1]
  mod_data <- mod_data[,-1]
  
  # Remove incomplete cases
  mod_data <- mod_data[complete.cases(mod_data), ]
  
  # Store column names
  CN <- names(mod_data)
  
  # Set default max value if not provided
  if (is.null(max)) {
    max <- floor(nrow(mod_data) / 5)
    cat("Maximum number of variables set to:", max, "\n")
  }
  
  # Scale the data (all columns except the last one, which is assumed to be the response)
  tryCatch({
    mod_data <- data.frame(cbind(scale(mod_data[,1:(ncol(mod_data) - 1)], TRUE, TRUE), 
                                 mod_data[[ncol(mod_data)]]))
    names(mod_data) <- CN
    row.names(mod_data) <- RN
  }, error = function(e) {
    stop(paste("Error during data scaling:", e$message))
  })
  
  # Separate leave-out samples if specified
  pred.data <- mod_data[row.names(mod_data) %in% leave.out, ]
  mod_data <- mod_data[!(row.names(mod_data) %in% leave.out), ]
  
  # Set cutoff for model selection
  cutoff <- 0.85
  
  # Run model search with correlation threshold
  cat("\nRunning model search with correlation threshold:", cor.threshold, "\n")
  tryCatch({
    models <- model.subset(data = mod_data, 
                           leave.out = leave.out, 
                           out.col = ncol(mod_data),
                           min = min, 
                           max = max, 
                           iterations = 1, 
                           cutoff = cutoff, 
                           cor.threshold = cor.threshold)
    
    # Print model table
    cat("\nTop models based on R-squared and Q-squared:\n")
    print(knitr::kable(models))
  }, error = function(e) {
    stop(paste("Error in model selection:", e$message))
  })
  
  # Let user select model if not provided
  if (is.null(what.model)) {
    what.model <- readline('Choose the model you would like to analyze (line number): ')
  }
  if (is.character(what.model)) {
    what.model <- as.numeric(what.model)
  }
  
  # Check if what.model is within valid range
  if (what.model < 1 || what.model > nrow(models)) {
    stop(paste("Invalid model selection. Please choose a number between 1 and", nrow(models)))
  }
  
  # Extract and display model coefficients
  tryCatch({
    mod.sum <- summary(lm(models$formula[what.model], mod_data))$coefficients
    cat('\nModel Coefficients:\n')
    colnames(mod.sum)[4] <- 'p value'
    print(knitr::kable(mod.sum))
  }, error = function(e) {
    stop(paste("Error extracting model coefficients:", e$message))
  })
  
  # Perform 3-fold cross-validation
  tryCatch({
    cv_3fold <- model.cv(models[what.model,1], mod_data, ncol(mod_data), 3, 50)
    dt3 <- data.frame(cv_3fold[[2]], cv_3fold[[1]])
    names(dt3) <- c('Q2', 'MAE')
    cat('\n3-fold CV Results:\n')
    print(knitr::kable(dt3))
  }, error = function(e) {
    warning(paste("Error in 3-fold cross-validation:", e$message))
    dt3 <- data.frame(Q2 = NA, MAE = NA)
  })
  
  # Perform 5-fold cross-validation
  tryCatch({
    cv_5fold <- model.cv(models[what.model,1], mod_data, ncol(mod_data), 5, 50)
    dt5 <- data.frame(cv_5fold[[2]], cv_5fold[[1]])
    names(dt5) <- c('Q2', 'MAE')
    cat('\n5-fold CV Results:\n')
    print(knitr::kable(dt5))
  }, error = function(e) {
    warning(paste("Error in 5-fold cross-validation:", e$message))
    dt5 <- data.frame(Q2 = NA, MAE = NA)
  })
  
  # Perform leave-one-out cross-validation
  tryCatch({
    cv_loo <- model.cv(models[what.model,1], mod_data, ncol(mod_data), nrow(mod_data), 1)
    dtloo <- data.frame(cv_loo[[2]], cv_loo[[1]])
    names(dtloo) <- c('Q2', 'MAE')
    cat('\nLOO-CV Results:\n')
    print(knitr::kable(dtloo))
  }, error = function(e) {
    warning(paste("Error in LOO cross-validation:", e$message))
    dtloo <- data.frame(Q2 = NA, MAE = NA)
  })
  
  # Make predictions for left-out samples if requested
  if (predict && nrow(pred.data) > 0) {
    tryCatch({
      model_formula <- models$formula[what.model]
      final_model <- lm(model_formula, mod_data)
      prediction <- predict(final_model, pred.data)
      real <- pred.data[[ncol(pred.data)]]
      
      # Calculate errors
      errors <- real - prediction
      prd.tab <- data.frame(
        "OOS Pred" = prediction, 
        "OOS Measured" = real, 
        "OOS Error" = errors
      )
      
      cat('\nOut-of-Sample Prediction Results:\n')
      print(knitr::kable(prd.tab))
      
      # Calculate and display RMSE
      rmse <- sqrt(mean(errors^2))
      mae <- mean(abs(errors))
      cat('\nOut-of-Sample Error Metrics:\n')
      print(knitr::kable(data.frame(RMSE = rmse, MAE = mae)))
      
      # Save predictions if requested
      if (save.pred) {
        pred_file <- paste0(dataset_name, "_model", what.model, "_predictions.csv")
        utils::write.csv(prd.tab, pred_file)
        cat("\nOut-of-sample predictions saved to:", pred_file, "\n")
      }
    }, error = function(e) {
      warning(paste("Error making out-of-sample predictions:", e$message))
      prd.tab <- NULL
    })
  }
  
  # Also show unnormalized coefficients
  tryCatch({
    mod_data_unn <- data.frame(data.table::fread(dataset, header = TRUE))
    mod.sum.unnormalized <- summary(lm(models$formula[what.model], mod_data_unn))$coefficients
    cat('\nUnnormalized Data Model Coefficients:\n')
    colnames(mod.sum.unnormalized)[4] <- 'p value'
    print(knitr::kable(mod.sum.unnormalized))
  }, error = function(e) {
    warning(paste("Error calculating unnormalized coefficients:", e$message))
  })
  
  # Prepare model visualization
  tryCatch({
    # Create table of statistics for plot annotation
    info.table <- data.frame(matrix(ncol = 1, nrow = 4))
    info.table[1,1] <- as.character(round(models$R.sq[what.model], 2))
    info.table[2,1] <- as.character(round(models$Q.sq[what.model], 2))
    info.table[3,1] <- as.character(round(dt5$Q2[1], 2))
    info.table[4,1] <- as.character(round(dt3$Q2[1], 2))
    row.names(info.table) <- c('R2', 'Q2_loo', 'Q2_5fold', 'Q2_3fold')
    names(info.table) <- 'stats'
    
    # Create text annotations for plot
    text1 <- paste(row.names(info.table)[1], info.table[1,1], sep = ' = ')
    text2 <- paste(row.names(info.table)[2], info.table[2,1], sep = ' = ')
    text3 <- paste(row.names(info.table)[3], info.table[3,1], sep = ' = ')
    text4 <- paste(row.names(info.table)[4], info.table[4,1], sep = ' = ')
    annotations <- stringr::str_c(c(text1, text2, text3, text4), collapse = "\n")
    
    # Create equation string for plot title
    equation <- list()
    equation[1] <- as.character(round(mod.sum[1,1], 2))  # Intercept
    for (i in 2:nrow(mod.sum)) {
      # Format each coefficient with its term
      if (mod.sum[i,1] < 0) {
        equation[i] <- paste(round(mod.sum[i,1], 2), row.names(mod.sum)[i], sep = ' * ')
      } else {
        equation[i] <- paste(paste('+', round(mod.sum[i,1], 2), sep = ''), 
                             row.names(mod.sum)[i], sep = ' * ')
      }
    }
    equation_str <- paste(equation, collapse = ' ')
    
    # Fit model and get prediction intervals
    best.mod <- lm(models$formula[what.model], mod_data)
    
    # Combine training and test data for plotting if prediction is requested
    if (predict && nrow(pred.data) > 0) {
      combined_data <- rbind(mod_data, pred.data)
      pred_interval <- predict(best.mod, newdata = combined_data, interval = 'pre', level = 0.9)
      plot.dat <- data.frame(
        "Measured" = combined_data[[ncol(combined_data)]], 
        "Predicted" = pred_interval[,1],
        "lwr" = pred_interval[,2],
        "upr" = pred_interval[,3]
      )
      rownames(plot.dat) <- c(rownames(mod_data), rownames(pred.data))
      
      # Mark which points are external test set
      plot.dat$is_test <- c(rep(FALSE, nrow(mod_data)), rep(TRUE, nrow(pred.data)))
    } else {
      # Only use training data
      pred_interval <- predict(best.mod, newdata = mod_data, interval = 'pre', level = 0.9)
      plot.dat <- data.frame(
        "Measured" = mod_data[[ncol(mod_data)]],
        "Predicted" = pred_interval[,1],
        "lwr" = pred_interval[,2],
        "upr" = pred_interval[,3]
      )
      rownames(plot.dat) <- rownames(mod_data)
      plot.dat$is_test <- FALSE
    }
    
    # Standardize row names formatting for chemical compounds
    row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat), "o_", '2-')
    row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat), "m_", '3-')
    row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat), "p_", '4-')
    row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat), "o4-", '2,4-')
    row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat), "m3-", '3,5-')
    row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat), "o3-", '2,3-')
    row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat), "basic", 'Ph')
    row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat), "penta_F", 'C6F5')
    
    # Add position classification based on common chemical naming patterns
    plot.dat$Position <- NA
    for (i in 1:nrow(plot.dat)) {
      if (plot.dat$is_test[i]) {
        plot.dat$Position[i] <- 'external'
      } else if (grepl('3-', row.names(plot.dat)[i])) {
        plot.dat$Position[i] <- 'meta'
      } else if (grepl('5-', row.names(plot.dat)[i])) {
        plot.dat$Position[i] <- 'meta'
      } else if (grepl('2-', row.names(plot.dat)[i])) {
        plot.dat$Position[i] <- 'ortho'
      } else if (grepl('Ph', row.names(plot.dat)[i])) {
        plot.dat$Position[i] <- 'Ph'
      } else if (grepl('C6F5', row.names(plot.dat)[i])) {
        plot.dat$Position[i] <- 'C6F5'
      } else if (grepl('4-', row.names(plot.dat)[i])) {
        plot.dat$Position[i] <- 'para'
      }
    }
    
    # Add labels and shape identifiers
    plot.dat$label <- row.names(plot.dat)
    plot.dat$shapes <- ifelse(plot.dat$is_test, 19, 18)
    
    # Calculate absolute errors for identifying outliers
    plot.dat$error <- abs(plot.dat$Measured - plot.dat$Predicted)
    
    # Model formula and dataset name for subtitles
    model_formula_str <- as.character(models$formula[what.model])
    plot_subtitle1 <- paste("Model:", model_formula_str)
    plot_subtitle2 <- paste("Dataset:", dataset_name)
    
    # Calculate plot limits
    x_range <- range(plot.dat$Measured, na.rm = TRUE)
    y_range <- range(plot.dat$Predicted, na.rm = TRUE)
    overall_range <- range(c(x_range, y_range), na.rm = TRUE)
    
    # Set error threshold for labeling outliers
    error_threshold <- 13  # Points with errors greater than this will be labeled
    
    # Create plot with ggplot2
    plot <- ggplot2::ggplot(plot.dat, ggplot2::aes(x = Measured, y = Predicted)) +
      # Points with color by position and shape by training/test
      ggplot2::geom_point(ggplot2::aes(color = Position, shape = as.factor(shapes)), size = 2) +
      
      # Prediction interval lines
      ggplot2::stat_smooth(ggplot2::aes(y = lwr), color = "cadetblue", linetype = "dashed",
                           se = FALSE, method = 'lm', fullrange = TRUE, linewidth = 0.8) +
      ggplot2::stat_smooth(ggplot2::aes(y = upr), color = "cadetblue", linetype = "dashed",
                           se = FALSE, method = 'lm', fullrange = TRUE, linewidth = 0.8) +
      
      # Regression line
      ggplot2::stat_smooth(method = 'lm', se = FALSE, formula = y ~ x,
                           color = 'black', fullrange = TRUE, linetype = 'dashed', linewidth = 0.8) +
      
      # Labels
      ggplot2::labs(
        x = 'Measured',
        y = 'Predicted',
        title = plot_title,
        subtitle = paste(plot_subtitle1, "\n", plot_subtitle2)
      ) +
      
      # Theme customization
      ggplot2::theme(
        axis.line.x = ggplot2::element_line(linewidth = 1, colour = "black"),
        axis.line.y = ggplot2::element_line(linewidth = 1, colour = "black"),
        axis.text.x = ggplot2::element_text(colour = "black", size = 12, face = 'bold'),
        axis.text.y = ggplot2::element_text(colour = "black", size = 12, face = 'bold'),
        axis.title.x = ggplot2::element_text(colour = "black", size = 12, face = 'bold'),
        axis.title.y = ggplot2::element_text(colour = "black", size = 12, face = 'bold'),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        
        # Legend positioning and styling
        legend.position = c(0.85, 0.15),
        legend.justification = c(0.5, 0.5),
        legend.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        legend.key.size = ggplot2::unit(0.8, "lines"),
        legend.title = ggplot2::element_text(size = 10),
        legend.text = ggplot2::element_text(size = 9),
        
        # Title and subtitle styling
        plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
        
        # Font family
        text = ggplot2::element_text(family = 'Arial')
      ) +
      
      # Color scheme for points
      ggplot2::scale_color_manual(
        '', 
        values = c(
          Ph = "black", 
          meta = 'tan1',
          para = '#66a182',
          ortho = '#d1495b', 
          external = 'steelblue4',
          C6F5 = 'darkgrey'
        )
      ) +
      
      # Shape mapping
      ggplot2::scale_shape_manual(
        values = c('18' = 18, '19' = 19), 
        guide = 'none'  # Hide the shape legend
      ) +
      
      # Set axis limits and aspect ratio
      ggplot2::coord_fixed(ratio = 1) +
      
      # Add text labels for outliers
      ggrepel::geom_text_repel(
        data = subset(plot.dat, error > error_threshold),
        ggplot2::aes(label = label),
        size = 3,
        min.segment.length = 0.5,
        seed = 42,
        point.padding = 0.4,
        segment.color = 'grey50',
        force_pull = 0.02,
        nudge_x = 0.03,
        direction = 'y'
      ) +
      
      # Add statistics annotation
      ggplot2::annotate(
        'text',
        x = min(plot.dat$lwr, na.rm = TRUE),
        y = max(plot.dat$Predicted, na.rm = TRUE), 
        label = annotations,
        parse = FALSE,
        hjust = "left", 
        vjust = 0
      )
    
    # If we have out-of-sample predictions, add a validation table to the plot
    if (predict && exists("prd.tab") && !is.null(prd.tab) && nrow(prd.tab) > 0) {
      prd.tab <- round(prd.tab, 2)
      # Create the validation table
      validation_table <- oos_validation_table(
        prd.tab,
        plot.title = paste("External Validation for", dataset_name),
        subtitle = paste("Model:", what.model),
        error_threshold = 20
      )
      
      # Print both plots
      gridExtra::grid.arrange(plot, validation_table, ncol = 2)
    } else {
      # Just return the main plot
      print(plot)
    }
    
    # Return the plot object
    invisible(plot)
    
  }, error = function(e) {
    warning(paste("Error creating plot:", e$message))
    # Return a minimal plot if there was an error
    ggplot2::ggplot() + 
      ggplot2::ggtitle("Error creating plot") +
      ggplot2::theme_minimal()
  })
}

#' Generate a Comprehensive Model Report with Visualizations from a Model List
#'
#' This function produces a detailed report for a linear model selected from a model list,
#' including cross-validation results, coefficient summaries, and visualizations. It handles
#' both training data evaluation and out-of-sample predictions, displaying results in both
#' tabular and graphical formats.
#'
#' @param dataset Path to a CSV file containing the data. First column should contain row identifiers,
#'        and the specified outcome column should contain the target variable.
#' @param model.list Path to a CSV file containing models generated by model.subset.parallel or similar
#'        functions. The file should include a 'formula' column with model formulas.
#' @param out.col Name of the outcome/target column in the dataset (default: 'output').
#' @param leave.out Vector of row identifiers to exclude from model training. These samples
#'        will be used for out-of-sample validation.
#' @param what.model Integer indicating which model from the model list to use (default: 1, the first model).
#' @param save.pred Logical indicating whether to save out-of-sample predictions to a CSV file (default: TRUE).
#' @param plot_title Title for the generated plot (default: "Linear Regression Model Analysis").
#'
#' @return A ggplot object visualizing the model performance, with both training and out-of-sample
#'         data points. The function also prints various statistics and tables to the console.
#'
#' @details
#' The function performs the following steps:
#' 1. Reads and preprocesses the dataset, scaling predictors and separating leave-out samples
#' 2. Extracts the specified model from the model list and evaluates it
#' 3. Performs cross-validation with different fold configurations (3-fold, 5-fold, and LOO)
#' 4. Creates a visualization of model performance showing predicted vs. measured values
#' 5. Makes predictions for left-out samples and calculates error metrics
#' 6. Generates both a standard scatter plot and a formatted validation table for out-of-sample results
#'
#' Special features of the visualization include:
#' - Color-coding points by position (ortho, meta, para, etc.)
#' - Displaying 90% prediction intervals
#' - Automated labeling of outlier points
#' - Summary statistics annotations (R², Q² for different CV schemes)
#'
#' @importFrom data.table fread
#' @importFrom knitr kable
#' @importFrom stringr str_c str_replace
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_point stat_smooth labs theme element_line element_text
#' @importFrom ggplot2 element_blank scale_color_manual xlim ylim coord_fixed annotate
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats lm predict complete.cases
#' @importFrom utils write.csv
#' @importFrom default default
#'
#' @export
model.report.from.list <- function(dataset,
                                   model.list,
                                   out.col = 'output',
                                   leave.out,
                                   what.model = 1,
                                   save.pred = T,
                                   plot_title = "Linear Regression Model Analysis") {
  default::default(data.frame) <- list(check.names = FALSE)
  mod_data <- data.frame(data.table::fread(dataset, header = T))
  RN <- mod_data[,1]
  mod_data <- mod_data[,-1]
  mod_data <- mod_data[complete.cases(mod_data), ]
  CN <- names(mod_data)
  out.col <- which(CN == out.col)
  mod_data <- data.frame(cbind(scale(mod_data[,-out.col], T, T), mod_data[, out.col]))
  names(mod_data)[1:(ncol(mod_data) - 1)] <- CN[-out.col]
  names(mod_data)[ncol(mod_data)] <- CN[out.col]
  row.names(mod_data) <- RN
  pred_data <- mod_data[row.names(mod_data) %in% leave.out, ]
  mod_data <- mod_data[!(row.names(mod_data) %in% leave.out), ]
  models <- data.frame(data.table::fread(model.list))
  mod.sum <- summary(lm(models$formula[what.model], mod_data))$coefficients
  cat('
  Model Coefficients')
  colnames(mod.sum)[4] <- 'p value'
  k.mod <- knitr::kable(mod.sum)
  print(k.mod)
  cv_3fold <- model.cv(models$formula[what.model], mod_data, dim(mod_data)[2], 3, 50)
  dt3 <- data.frame(cv_3fold[[2]], cv_3fold[[1]])
  names(dt3) <- c('Q2', 'MAE')
  cat('
  3-fold CV')
  tab_dt3 <- knitr::kable(dt3)
  print(tab_dt3)
  
  cv_5fold <- model.cv(models$formula[what.model], mod_data, dim(mod_data)[2], 5, 50)
  dt5 <- data.frame(cv_5fold[[2]], cv_5fold[[1]])
  names(dt5) <- c('Q2', 'MAE')
  cat('
  5-fold CV')
  tab_dt5 <- knitr::kable(dt5)
  print(tab_dt5)
  
  cv_loo <- model.cv(models$formula[what.model], mod_data, dim(mod_data)[2], nrow(mod_data), 1)
  dtloo <- data.frame(cv_loo[[2]], cv_loo[[1]])
  names(dt3) <- c('Q2', 'MAE')
  cat('
  LOO-CV')
  tab_dtloo <- knitr::kable(dtloo)
  print(tab_dtloo)
  
  
  mod_data_unn <- data.frame(data.table::fread(dataset, header = T))
  mod.sum.unnormalized <- summary(lm(models$formula[what.model], mod_data_unn))$coefficients
  cat('
  Unnormalized Data Model Coefficients')
  colnames(mod.sum.unnormalized)[4] <- 'p value'
  k.mod.unn <- knitr::kable(mod.sum.unnormalized)
  print(k.mod.unn)
  
  ## model.plot
  info.table <- data.frame(matrix(ncol = 1, nrow = 4))
  info.table[1,1] <- as.character(round(models$R.sq[what.model], 2))
  info.table[2,1] <- as.character(round(dtloo[1, 1], 2))
  info.table[3,1] <- as.character(round(dt5[1, 1], 2))
  info.table[4,1] <- as.character(round(dt3[1, 1], 2))
  row.names(info.table) <-  c('R2', 'Q2_loo', 'Q2_5fold', 'Q2_3fold')
  names(info.table) <- 'stats'
  text1 <- paste(row.names(info.table)[1], info.table[1,1], sep = ' = ')
  text2 <- paste(row.names(info.table)[2], info.table[2,1], sep = ' = ')
  text3 <- paste(row.names(info.table)[3], info.table[3,1], sep = ' = ')
  text4 <- paste(row.names(info.table)[4], info.table[4,1], sep = ' = ')
  
  annotations <- stringr::str_c(c(text1,
                                  text2,
                                  text3,
                                  text4),
                                collapse = "\n")
  
  model = models$formula[what.model]
  best.mod <- lm(model, data = mod_data)
  pred_interval <- predict(best.mod,
                           newdata = mod_data,
                           interval = 'pre',
                           level = 0.9)
  plot.dat <- data.frame(cbind(mod_data[dim(mod_data)[2]], pred_interval))
  colnames(plot.dat) <- c('Measured', 'Predicted', 'lwr', 'upr')
  rownames(plot.dat) <- row.names(mod_data)
  
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o_",'2-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"m_",'3-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"p_",'4-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o4-",'2,4-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"m3-",'3,5-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o3-",'2,3-')
  
  plot.dat <- dplyr::mutate(plot.dat, Position = rep(NA, nrow(plot.dat)))
  
  for (i in 1:nrow(mod_data)) {
    if (grepl('3-',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'meta'
    }
    if (grepl('5-',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'meta'
    }
    if (grepl('2-',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'ortho'
    }
    if (grepl('basic',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'Ph'
    }
    if (grepl('penta_F',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'C6F5'
    }
    if (grepl('4-',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'para'
    }
  }
  plot.dat <- dplyr::mutate(plot.dat, label = row.names(plot.dat))
  
  plot.dat <- dplyr::mutate(plot.dat, 
                            shapes = c(rep(18, nrow(mod_data))))
  
  
  # Calculate absolute errors
  plot.dat$error <- abs(plot.dat$Measured - plot.dat$Predicted)
  # Create label vector that only includes points with large errors
  plot.dat$show_label <- ifelse(plot.dat$error > 15, plot.dat$label, "")
  
  # Extract dataset name (without path and extension) for subtitle
  dataset_name <- basename(dataset)
  dataset_name <- sub("\\.[^.]*$", "", dataset_name)
  
  # Get model formula as string for subtitle
  model_formula_str <- as.character(models$formula[what.model])
  
  # Store for plot subtitles
  plot_subtitle1 <- paste("Model:", model_formula_str)
  plot_subtitle2 <- paste("Dataset:", dataset_name)
  
  # Calculate plot limits for proper legend positioning and for x/y limits
  x_min <- min(plot.dat[1:nrow(mod_data),3])
  x_max <- max(plot.dat[1:nrow(mod_data),4])
  y_min <- min(plot.dat[1:nrow(mod_data),3])
  y_max <- max(plot.dat[1:nrow(mod_data),4])
  
  plot <- suppressMessages(ggplot2::ggplot(plot.dat, ggplot2::aes(x = Measured, y = Predicted)) +
                             ggplot2::geom_point(size = 2, shape = plot.dat$shapes, ggplot2::aes(color = Position)) +
                             ggplot2::stat_smooth(ggplot2::aes(y = lwr), color = "cadetblue", linetype = "dashed",
                                                  se = F, method = 'lm', fullrange = T, size = 0.8) +
                             ggplot2::stat_smooth(ggplot2::aes(y = upr), color = "cadetblue", linetype = "dashed",
                                                  se = F, method = 'lm', fullrange = T, size = 0.8) +
                             ggplot2::labs(x = 'Measured',
                                           y = 'Predicted',
                                           title = plot_title, 
                                           subtitle = paste(plot_subtitle1, "\n", plot_subtitle2)) +
                             ggplot2::stat_smooth(method = 'lm',se = F, formula = y~x,
                                                  color = 'black',fullrange = T, linetype = 'dashed') +
                             ggplot2::theme(axis.line.x = ggplot2::element_line(linewidth = 1, colour = "black"),
                                            axis.line.y = ggplot2::element_line(linewidth = 1, colour = "black"),
                                            axis.text.x = ggplot2::element_text(colour = "black", size = 12,face = 'bold'),
                                            axis.text.y = ggplot2::element_text(colour = "black", size = 12,face = 'bold'),
                                            axis.title.x = ggplot2::element_text(colour = "black", size = 12,face = 'bold'),
                                            axis.title.y = ggplot2::element_text(colour = "black", size = 12,face = 'bold'),
                                            panel.grid.major = ggplot2::element_blank(),
                                            panel.grid.minor = ggplot2::element_blank(),
                                            panel.border = ggplot2::element_blank(), 
                                            panel.background = ggplot2::element_blank(),
                                            # Position legend at the bottom right inside the plot
                                            legend.position = c(0.85, 0.15),
                                            legend.justification = c(0.5, 0.5),
                                            legend.background = ggplot2::element_blank(),
                                            legend.key = ggplot2::element_blank(),
                                            legend.key.size = unit(0.8, "lines"),
                                            legend.title = ggplot2::element_text(size = 10),
                                            legend.text = ggplot2::element_text(size = 9),
                                            # Title and subtitle styling
                                            plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                                            plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5)) +
                             ggplot2::scale_color_manual('', values = c(Ph = "black", meta = 'tan1', C6F5 = 'darkgrey',
                                                                        para = '#66a182',ortho = '#d1495b', external = 'steelblue4')) +
                             ggplot2::xlim(x_min, x_max) +
                             ggplot2::ylim(y_min, y_max) +
                             ggplot2::coord_fixed(ratio = 1) +
                             ggrepel::geom_text_repel(
                               data = subset(plot.dat, error > 13),
                               aes(label = label),
                               size = 3,
                               min.segment.length = Inf,
                               seed = 42,
                               point.padding = 0.4,
                               segment.color = 'grey50',
                               force_pull = 0.02,
                               nudge_x = 0.022,
                               direction = 'y'
                             ) +
                             ggplot2::theme(text = ggplot2::element_text(family = 'Arial')) +
                             ggplot2::annotate('text',
                                               x = min(plot.dat[1:nrow(mod_data),3]),
                                               y = max(plot.dat[1:nrow(mod_data),2]), label = annotations,
                                               parse = F,
                                               hjust = "left", vjust = 0))
  plot
  
  
  prediction <- round(predict(lm(models$formula[what.model], mod_data), pred_data), 0)
  real <- pred_data[, dim(mod_data)[2]]
  RMSE <- round(sqrt((real-prediction)^2), 0)
  prd.tab <- data.frame(prediction, real, RMSE)
  names(prd.tab) <- c('OOS Pred', 'OOS Measured', 'OOS Error')
  if (save.pred == T) write.csv(prd.tab, "OOS_predictions.csv", row.names = TRUE)
  k.prd.tab <- knitr::kable(prd.tab)
  print(k.prd.tab)
  
  print(knitr::kable(data.frame('Mean of RMSE' = mean(RMSE))))
  
  oos_data <- prd.tab[, 1:2]
  
  row.names(oos_data) <- stringr::str_replace(row.names(oos_data),"o_",'2-')
  row.names(oos_data) <- stringr::str_replace(row.names(oos_data),"m_",'3-')
  row.names(oos_data) <- stringr::str_replace(row.names(oos_data),"p_",'4-')
  row.names(oos_data) <- stringr::str_replace(row.names(oos_data),"o4-",'2,4-')
  row.names(oos_data) <- stringr::str_replace(row.names(oos_data),"m3-",'3,5-')
  row.names(oos_data) <- stringr::str_replace(row.names(oos_data),"o3-",'2,3-')
  row.names(oos_data) <- stringr::str_replace(row.names(oos_data),"basic",'Ph')
  row.names(oos_data) <- stringr::str_replace(row.names(oos_data),"penta_F",'C6F5')
  
  print(oos_validation_table(prd.tab))
  
  # Final plot with out-of-sample points
  plot + 
    ggplot2::geom_point(data = oos_data,
                        ggplot2::aes(x = `OOS Measured`, y = `OOS Pred`,
                                     color = "Out-of-Sample"),  # Added color aesthetic
                        size = 1.7,
                        shape = 19) +
    ggplot2::scale_color_manual('', values = c(Ph = "black", 
                                               C6F5 = "darkgrey",
                                               meta = 'tan1',
                                               para = '#66a182', 
                                               ortho = '#d1495b', 
                                               external = 'steelblue4',
                                               "Out-of-Sample" = "cadetblue3")) +  # Added out-of-sample color
    ggrepel::geom_text_repel(data = oos_data,
                             ggplot2::aes(x = `OOS Measured`, y = `OOS Pred`, 
                                          label = rownames(oos_data)),
                             size = 3,
                             min.segment.length = Inf,
                             seed = 42,
                             point.padding = 0.4,
                             segment.color = 'cadetblue3',
                             force_pull = 0.01,
                             nudge_y = -0.022,
                             direction = 'y')
}


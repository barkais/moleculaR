####### ----------------------------------------------------#####
####### -----------------Utility Functions------------------#####
####### ----------------------------------------------------#####

#' Cross validate (k-fold) a single model using parallel computation
#'
#' An iterative version of model.single.cv
#' @param formula model formula (not designed for direct use)
#' @param data data frame with outcome column
#' @param out.col column number with output
#' @param folds define K
#' @param iterations Number of iterations
#' @export
#' @return averaged MAE and Q2
model.cv.parallel <- function(formula, data, out.col, folds, iterations) {
  model.single.cv.iterator <- function(i) {
    tool <- model.single.cv(formula, data, out.col, folds)
    return(data.frame(tool[[1]], tool[[2]]))
  }
  results <- do.call(rbind, parallel::mclapply(1:iterations, model.single.cv.iterator))
  MAE.validation <- Reduce(`+`, results[, 1]) / iterations
  q2.validation <- Reduce(`+`, results[, 2]) / iterations
  return(list(MAE.validation, q2.validation))
}

#' Parallelized Subset Feature Selection with Correlation Threshold
#'
#' This function generates all possible subsets of features in a user-defined range (min to max),
#' ensuring that no two variables in a subset have a Pearson correlation coefficient above the
#' specified correlation threshold (`cor.threshold`). It performs parallel computation for faster execution.
#'
#' @param data A data frame containing the features and the target column.
#' @param out.col An integer indicating the column index of the target variable (default: last column).
#' @param min The minimum number of features to include in a subset (default: 2).
#' @param max The maximum number of features to include in a subset (default: floor(nrow(data) / 5)).
#' @param results_name Name of file containing results. The function writes results on the fly, to free up RAM.
#' @param folds Number of folds for cross-validation (default: number of rows in `data`).
#' @param iterations Number of iterations for cross-validation (default: 1).
#' @param cutoff R-squared cutoff for model selection (default: 0.85).
#' @param cor.threshold Pearson correlation threshold for feature selection (default: 0.7).
#' @return A data frame with the best models, their R-squared, Q-squared, and MAE values.
#' @export

model.subset.parallel <- function(data, out.col = dim(data)[2],
                                  min = 2, max = floor(dim(data)[1] / 5),
                                  results_name = 'results_df',
                                  folds = nrow(data), iterations = 1,
                                  cutoff = 0.85, cor.threshold = 0.7) {
  
  # Helper function to process results files in batches
  process_results_batch <- function(file_list, batch_size=50) {
    best_models <- data.frame()
    
    for(i in seq(1, length(file_list), by=batch_size)) {
      end_idx <- min(i + batch_size - 1, length(file_list))
      current_files <- file_list[i:end_idx]
      
      # Read and process batch
      batch_results <- lapply(current_files, function(f) {
        res <- try({
          data <- data.table::fread(f)
          file.remove(f)  # Clean up immediately
          return(data)
        }, silent = TRUE)
        
        if(inherits(res, "try-error")) {
          warning(paste("Error reading file:", f))
          return(NULL)
        }
        return(res)
      })
      
      # Remove any NULL results from failed reads
      batch_results <- Filter(Negate(is.null), batch_results)
      
      if(length(batch_results) > 0) {
        # Combine and filter
        batch_df <- do.call(rbind, batch_results)
        batch_df <- dplyr::arrange(batch_df, desc(V2))  # V2 is R.sq
        if(nrow(batch_df) > 50) batch_df <- batch_df[1:50,]
        
        # Update best models
        best_models <- rbind(best_models, batch_df)
        best_models <- dplyr::arrange(best_models, desc(V2))
        if(nrow(best_models) > 50) best_models <- best_models[1:50,]
      }
      
      # Clear memory
      rm(batch_results, batch_df)
      gc()
    }
    
    return(best_models)
  }
  
  # Get the output variable and feature names
  output <- stringr::str_c("`", names(data[out.col]), "`")
  vars <- names(data[, -out.col])
  for (i in 1:length(vars)) {
    vars[i] <- stringr::str_c("`", vars[i], "`")
  }
  
  # Compute the correlation matrix
  cor_matrix <- cor(data[, -out.col])
  
  # Process each feature size in range
  for (i in min:max) {
    # Calculate combinations in chunks to save memory
    chunk_size <- min(1000000, choose(length(vars), i))  # Adjust chunk_size based on available memory
    n_vars <- length(vars)
    total_combs <- choose(n_vars, i)
    
    for(chunk_start in seq(1, total_combs, by=chunk_size)) {
      chunk_end <- min(chunk_start + chunk_size - 1, total_combs)
      
      # Generate only this chunk of combinations
      current_chunk <- combn(vars, i, simplify=FALSE)[chunk_start:chunk_end]
      
      # Process combinations in this chunk
      chunk_df <- data.frame(
        formula = sapply(current_chunk, function(x) paste(x, collapse=" + ")),
        stringsAsFactors = FALSE
      )
      
      # Apply correlation filtering
      filtered_chunk <- sapply(chunk_df$formula, function(formula) {
        features <- strsplit(formula, " + ")[[1]]
        feature_indices <- which(colnames(cor_matrix) %in% gsub("`", "", features))
        if(length(feature_indices) > 1) {
          feature_cor <- cor_matrix[feature_indices, feature_indices]
          !any(abs(feature_cor[upper.tri(feature_cor)]) > cor.threshold)
        } else {
          TRUE
        }
      })
      
      chunk_df <- chunk_df[filtered_chunk,]
      
      if(nrow(chunk_df) > 0) {
        # Add output variable to formulas
        chunk_df$formula <- paste(output, "~", chunk_df$formula)
        
        # Calculate R-squared in parallel
        r_squared_results <- parallel::mclapply(
          chunk_df$formula,
          function(x) try(summary(lm(x, data = data))$r.squared, silent=TRUE),
          mc.cores = parallel::detectCores() - 1
        )
        
        # Filter out any errors
        valid_results <- !sapply(r_squared_results, inherits, "try-error")
        chunk_df <- chunk_df[valid_results,]
        r_squared_results <- unlist(r_squared_results[valid_results])
        
        if(length(r_squared_results) > 0) {
          chunk_df$R.sq <- r_squared_results
          
          # Keep only the best models from this chunk
          chunk_df <- chunk_df[chunk_df$R.sq > cutoff,]
          chunk_df <- dplyr::arrange(chunk_df, desc(R.sq))
          if(nrow(chunk_df) > 10) chunk_df <- chunk_df[1:10,]
          
          # Write results if we have any
          if(nrow(chunk_df) > 0) {
            write.table(chunk_df,
                        file = paste0(results_name, format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
                        append = TRUE, col.names = FALSE)
          }
        }
      }
      
      # Clear memory
      rm(current_chunk, chunk_df)
      gc()
    }
  }
  
  # Process all result files
  file_list <- list.files(pattern = paste0(results_name, "\\d{8}_\\d{6}\\.csv$"))
  
  if(length(file_list) == 0) {
    stop("No models met the criteria")
  }
  
  # Process results in batches
  forms.cut <- process_results_batch(file_list)
  names(forms.cut) <- c('formula', 'R.sq')
  
  # Select top 50 models for cross-validation
  if(nrow(forms.cut) > 50) forms.cut <- forms.cut[1:50,]
  
  # Perform cross-validation on the top models
  q2.list <- list()
  mae.list <- list()
  
  for(i in 1:nrow(forms.cut)) {
    stts <- model.cv.parallel(formula = forms.cut$formula[i],
                              data = data,
                              out.col = out.col,
                              folds = folds,
                              iterations = iterations)
    q2.list[[i]] <- stts[2]
    mae.list[[i]] <- stts[1]
  }
  
  # Add Q-squared and MAE to the final results
  forms.cut$Q.sq <- unlist(q2.list)
  forms.cut$MAE <- unlist(mae.list)
  
  # Sort by Q-squared and add model numbers
  forms.cut <- dplyr::arrange(forms.cut, desc(Q.sq))
  forms.cut$Model <- seq_len(nrow(forms.cut))
  
  # Return top 15 models
  return(forms.cut[1:15,])
}

####### ----------------------------------------------------#####
####### -------------------User Functions-------------------#####
####### ----------------------------------------------------#####

#'  Generate a .csv file with top 10 linear models resulting from a parallel 
#'  model.subset search. 
#'
#' Screen for models and cross validate them.
#' @param dataset a dataframe with outcome column (must be named 'output')
#' @param min minimum # of features (default = 2)
#' @param max max # of features (defaults = # of observations / 5)
#' @param results_name Name of file containing results. The function writes results on the fly, to free up RAM.
#' @param leave.out name of observations to leave out (e.g. 'p_Br')
#' @param folds number of folds for CV
#' @param iterations number of iterations for CV
#' @return models list, CV results for k=3/5/LOO.
#' @aliases models.list.parallel
#' @export
models.list.parallel <- function(dataset,
                                 min = 2,
                                 max = floor(dim(mod_data)[1] / 5),
                                 results_name = 'results_df',
                                 leave.out = '',
                                 folds = nrow(mod_data), 
                                 iterations = 1) {
  
  # Create directory name based on parameters and timestamp
  dir_name <- paste0(
    tools::file_path_sans_ext(basename(dataset)),
    '_min', as.character(min),
    '_max', as.character(max),
    '_', format(Sys.time(), "%Y%m%d_%H%M%S")
  )
  
  # Create directory if it doesn't exist
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  # Modify results_name to include directory path
  results_name <- file.path(dir_name, results_name)
  
  default::default(data.frame) <- list(check.names = FALSE)
  cat(tools::file_path_sans_ext(basename(dataset)))
  mod_data <- data.frame(data.table::fread(dataset, header = T),
                         check.names = F)
  RN <- mod_data[,1]
  mod_data <- mod_data[,-1]
  mod_data <- mod_data[complete.cases(mod_data), ]
  CN <- names(mod_data)
  mod_data <- data.frame(cbind(scale(mod_data[,1:dim(mod_data)[2] - 1], T, T), mod_data$output))
  names(mod_data) <- CN
  row.names(mod_data) <- RN
  pred.data <- mod_data[row.names(mod_data) %in% leave.out, ]
  mod_data <- mod_data[!(row.names(mod_data) %in% leave.out), ]
  
  models <- model.subset.parallel(mod_data,
                                  min = min,
                                  max = max,
                                  results_name = results_name,
                                  folds = folds, 
                                  iterations = iterations)
  
  tab <- knitr::kable(models)
  print(tab)
  
  # Save models list in the same directory
  output_file <- file.path(
    dir_name,
    paste0(
      tools::file_path_sans_ext(basename(dataset)),
      '_',
      as.character(min),
      '_',
      as.character(max),
      '_models.list.csv'
    )
  )
  
  write.csv(models, output_file)
  
  # Instead of deleting, just return the directory name
  cat("\nResults saved in directory:", dir_name, "\n")
  
  return(models)
}

#'  Generate models and a plot. A parallel computation version.
#'
#' Screen for models, cross validate them and plot. Designed for interactive work.
#' @param dataset a dataframe with outcome column (must be named 'output')
#' @param min minimum # of features (default = 2)
#' @param max max # of features (defaults = # of observations / 5)
#' @param leave.out name of observations to leave out (e.g. 'p_Br')
#' @param predict if leave.out is not empty, should a prediction
#' of it be computed?
#' @param ext.val if leave.out is not empty, should this list be used as a test?
#' @param what.model if not in an interactive session, what model should be used?
#' @return models list, CV results for k=3/5/LOO and a plot of choice. interactive.
#' @export
model.report.parallel <- function(dataset, min = 2, max = floor(dim(mod_data)[1] / 5),
                         leave.out = '', predict = F, ext.val = F,
                         what.model = NULL) {
  default::default(data.frame) <- list(check.names = FALSE)
  cat(tools::file_path_sans_ext(basename(dataset)))
  mod_data <- data.frame(data.table::fread(dataset, header = T))
  RN <- mod_data[,1]
  mod_data <- mod_data[,-1]
  mod_data <- mod_data[complete.cases(mod_data), ]
  CN <- names(mod_data)
  mod_data <- data.frame(cbind(scale(mod_data[,1:dim(mod_data)[2] - 1], T, T), mod_data$output))
  names(mod_data) <- CN
  row.names(mod_data) <- RN
  pred.data <- mod_data[row.names(mod_data) %in% leave.out, ]
  mod_data <- mod_data[!(row.names(mod_data) %in% leave.out), ]
  cutoff <- ifelse(ext.val == T, 0.8, 0.9)
  models <- model.subset.parallel(mod_data,
                                  min = min,
                                  max = max,
                                  iterations = 1,
                                  cutoff = cutoff)
  if (ext.val == T) {
    MAE.list <- vector(length = nrow(models))
    for (model in 1:nrow(models)) {
      prediction <- predict(lm(models[model, 1], mod_data), pred.data)
      real <- pred.data$output
      MAE.list[model] <- mean(abs(prediction - real))
    }
    models$MAE <- MAE.list
    models$Model <- rep(0, nrow(models))
    models <- dplyr::arrange(models, models$MAE)
    models$Model <- 1:nrow(models)
    names(models)[4] <- 'Pred.MAE'
  }
  tab <- knitr::kable(models)
  print(tab)
  if (is.null(what.model)) what.model <- readline('Choose the model you would like to plot (line number): ')
  if (is.character(what.model)) what.model <- as.numeric(what.model)
  mod.sum <- summary(lm(models[what.model, 1], mod_data))$coefficients
  cat('
  Model Coefficients')
  colnames(mod.sum)[4] <- 'p value'
  k.mod <- knitr::kable(mod.sum)
  print(k.mod)
  cv_3fold <- model.cv.parallel(models[what.model,1], mod_data, dim(mod_data)[2], 3, 50)
  dt3 <- data.frame(cv_3fold[[2]], cv_3fold[[1]])
  names(dt3) <- c('Q2', 'MAE')
  cat('
  3-fold CV')
  tab_dt3 <- knitr::kable(dt3)
  print(tab_dt3)
  cv_5fold <- model.cv.parallel(models[what.model,1], mod_data, dim(mod_data)[2], 5, 50)
  dt5 <- data.frame(cv_5fold[[2]], cv_5fold[[1]])
  names(dt5) <- c('Q2', 'MAE')
  cat('
  5-fold CV')
  tab_dt5 <- knitr::kable(dt5)
  print(tab_dt5)
  if (predict == T) {
    prediction <- predict(lm(models[what.model, 1], mod_data), pred.data)
    real <- pred.data$output
    prd.tab <- data.frame(prediction, real)
    names(prd.tab) <- c('OOS Pred', 'OOS Measured')
    k.prd.tab <- knitr::kable(prd.tab)
    print(k.prd.tab)
  }
  mod_data_unn <- data.frame(data.table::fread(dataset, header = T, check.names = T))
  mod.sum.unnormalized <- summary(lm(models[what.model, 1], mod_data_unn))$coefficients
  cat('
  Unnormalized Data Model Coefficients')
  colnames(mod.sum.unnormalized)[4] <- 'p value'
  k.mod.unn <- knitr::kable(mod.sum.unnormalized)
  print(k.mod.unn)
  
  ## model.plot
  info.table <- data.frame(matrix(ncol = 1, nrow = 4))
  info.table[1,1] <- as.character(round(models[what.model, 2], 2))
  info.table[2,1] <- as.character(round(models[what.model, 3], 2))
  info.table[3,1] <- as.character(round(dt5[1, 1], 2))
  info.table[4,1] <- as.character(round(dt3[1, 1], 2))
  row.names(info.table) <-  c('R2', 'Q2_loo', 'Q2_5fold', 'Q2_3fold')
  names(info.table) <- 'stats'
  text1 <- paste(row.names(info.table)[1], info.table[1,1], sep = ' = ')
  text2 <- paste(row.names(info.table)[2], info.table[2,1], sep = ' = ')
  text3 <- paste(row.names(info.table)[3], info.table[3,1], sep = ' = ')
  text4 <- paste(row.names(info.table)[4], info.table[4,1], sep = ' = ')
  
  equation <- list()
  for (i in 2:nrow(mod.sum)) {
    equation[i] <- ifelse(mod.sum[i,1] < 0, paste(round(mod.sum[i,1], 2),
                                                  row.names(mod.sum)[i], sep = ' * '),
                          paste(paste('+', round(mod.sum[i,1], 2), sep = ''),
                                row.names(mod.sum)[i], sep = ' * '))
  }
  equation[1] <- as.character(round(mod.sum[1,1], 2))
  equation <- paste(equation,  collapse = ' ')
  
  annotations <- stringr::str_c(c(text1,
                                  text2,
                                  text3,
                                  text4),
                                collapse = "\n")
  
  model = models[what.model, 1]
  data = data.frame(rbind(mod_data, pred.data))
  best.mod <- lm(model, data = mod_data)
  pred_interval <- predict(best.mod,
                           newdata = data,
                           interval = 'pre',
                           level = 0.9)
  plot.dat <- data.frame(cbind(data$output, pred_interval))
  colnames(plot.dat) <- c('Measured', 'Predicted', 'lwr', 'upr')
  rownames(plot.dat) <- row.names(data)
  lo.names <- which((row.names(plot.dat) %in% leave.out))
  
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o_",'2-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"m_",'3-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"p_",'4-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o4-",'2,4-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"m3-",'3,5-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o3-",'2,3-')
  
  plot.dat <- dplyr::mutate(plot.dat, Position = rep(NA, nrow(plot.dat)))
  plot.dat$Position[lo.names] <- 'external'
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
    if (grepl('4-',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'para'
    }
  }
  plot.dat <- dplyr::mutate(plot.dat, label = row.names(plot.dat))
  plot.dat <- dplyr::mutate(plot.dat, 
                            shapes = c(rep(18, nrow(mod_data)),
                                       rep(19, nrow(pred.data))))
  plot <- ggplot2::ggplot(plot.dat, ggplot2::aes(x = Measured, y = Predicted)) +
    ggplot2::geom_point(size = 2, shape = plot.dat$shapes, ggplot2::aes(color = Position)) +
    ggplot2::stat_smooth(ggplot2::aes(y = lwr), color = "cadetblue", linetype = "dashed",
                         se = F, method = 'lm', fullrange = T, size = 0.8) +
    ggplot2::stat_smooth(ggplot2::aes(y = upr), color = "cadetblue", linetype = "dashed",
                         se = F, method = 'lm', fullrange = T, size = 0.8) +
    ggplot2::labs(x = 'Measured',y = 'Predicted') +
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
                   panel.border = ggplot2::element_blank(), panel.background = ggplot2::element_blank(),
                   legend.position = c(2,2)) +
    ggplot2::scale_color_manual('', values = c(Ph = "black", meta = 'tan1',
                                               para = '#66a182',ortho = '#d1495b', external = 'steelblue4')) +
    ggplot2::xlim(min(plot.dat[,3]), max(plot.dat[,4])) +
    ggplot2::ylim(min(plot.dat[,3]), max(plot.dat[,4])) +
    ggplot2::coord_fixed(ratio = 1) +
    ggrepel::geom_text_repel(size = 3,
                             ggplot2::aes(label = label),
                             min.segment.length = Inf,
                             seed = 42,
                             point.padding = 0.4,
                             segment.color = 'white',
                             force_pull = 0.02,
                             nudge_x = 0.022,
                             direction = 'y') +
    ggplot2::theme(text = ggplot2::element_text(family = 'Arial'),
                   plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = -0.3)) +
    ggplot2::annotate('text',
                      x = min(plot.dat[,3]),
                      y = max(plot.dat[,2]), label = annotations,
                      parse = F,
                      hjust = "left", vjust = 0) +
    ggplot2::ggtitle(equation)
  plot
}


####### ----------------------------------------------------#####
####### -----------------Utility Functions------------------#####
####### ----------------------------------------------------#####

#' @keywords internal
#'
#' @description
#' OS detection and function swapping to provide seamless cross-platform compatibility.
#' On Windows systems, this code replaces the standard model.subset.parallel function
#' with the Windows-compatible version (model.subset.parallel.windows).
#'
#' This enables package users to always call model.subset.parallel() regardless of
#' their operating system, with the appropriate implementation being used automatically.
#'
#' @noRd
if(.Platform$OS.type == "windows") {
  # Store the original function for reference
  model.subset.parallel.original <- model.subset.parallel
  
  # Replace with Windows-compatible version
  model.subset.parallel <- model.subset.parallel.windows
}

#' Parallelized Subset Feature Selection with Correlation Threshold
#'
#' This function generates all possible subsets of features in a user-defined range (min to max),
#' ensuring that no two variables in a subset have a Pearson correlation coefficient above the
#' specified correlation threshold (`cor.threshold`). It performs parallel computation for faster execution
#' and handles large feature spaces efficiently by processing combinations in chunks and writing
#' intermediate results to temporary files.
#'
#' @param data A data frame containing the features and the target column.
#' @param out.col An integer indicating the column index of the target variable (default: last column).
#' @param min The minimum number of features to include in a subset (default: 2).
#' @param max The maximum number of features to include in a subset (default: floor(nrow(data) / 5)).
#' @param results_name Base name for temporary files storing model results (default: 'results_df').
#' @param folds Number of folds for cross-validation (default: number of rows in `data`).
#' @param iterations Number of iterations for cross-validation (default: 1).
#' @param cutoff Initial R-squared cutoff for model selection (default: 0.65). The function adaptively adjusts
#'        this value if too few models meet the initial threshold.
#' @param cor.threshold Pearson correlation threshold for feature selection (default: 1). Variable pairs with
#'        absolute correlation exceeding this value will not appear together in the same model.
#' @param min_models Minimum number of models to retain for cross-validation (default: 50). The function
#'        adaptively adjusts the cutoff to ensure at least this many models are evaluated.
#' @param verbose Logical indicating whether to print progress information (default: TRUE).
#'
#' @return A data frame with the best models (up to 15), ranked by cross-validated Q-squared value. 
#'         The data frame contains columns for each model's formula, R-squared, Q-squared, MAE, and 
#'         a Model number.
#'
#' @details
#' The function works in several stages:
#' 1. For each feature subset size from `min` to `max`, it generates combinations and filters them
#'    based on the correlation threshold.
#' 2. It evaluates models for these combinations in parallel, saving results to temporary files.
#' 3. If too few models meet the initial R-squared cutoff, it adaptively lowers the threshold.
#' 4. It performs cross-validation on the top models to compute Q-squared and MAE values.
#' 5. It ranks the final models by Q-squared and returns the top 15.
#'
#' The function handles memory constraints by processing combinations in chunks and writing
#' intermediate results to disk. This approach allows it to explore very large feature spaces
#' that would otherwise exhaust available RAM.
#'
#' @importFrom stringr str_c
#' @importFrom parallel mclapply detectCores
#' @importFrom utils write.csv read.csv combn
#' @importFrom stats cor lm
#'
#' @export

model.subset.parallel <- function(data, out.col = dim(data)[2],
                                  min = 2, max = floor(dim(data)[1] / 5),
                                  results_name = 'results_df',
                                  folds = nrow(data), iterations = 1,
                                  cutoff = 0.65, cor.threshold = 1,
                                  min_models = 50, verbose = TRUE) {
  
  # Get the output variable and feature names
  output <- stringr::str_c("`", names(data[out.col]), "`")
  vars <- names(data[, -out.col])
  for (i in 1:length(vars)) {
    vars[i] <- stringr::str_c("`", vars[i], "`")
  }
  
  # Compute the correlation matrix for features
  cor_matrix <- cor(data[, -out.col])
  
  # Extract directory path from results_name
  results_dir <- dirname(results_name)
  if(results_dir == ".") results_dir <- getwd()
  
  # Clean up any existing result files with the same name pattern in this directory
  file_pattern <- paste0('^', basename(results_name), ".*\\.csv$")
  old_files <- list.files(path = results_dir, pattern = file_pattern, full.names = TRUE)
  if(length(old_files) > 0) {
    if(verbose) cat("Cleaning up", length(old_files), "old result files...\n")
    for(file in old_files) {
      unlink(file)
    }
  }
  
  # Track total models evaluated
  total_models_evaluated <- 0
  models_saved <- 0
  
  # Process each feature size in range
  for (i in min:max) {
    if(verbose) cat("Processing feature size:", i, "\n")
    
    # Safety check for combinations
    if(choose(length(vars), i) <= 0 || i > length(vars)) {
      if(verbose) cat("  Skipping - invalid feature size\n")
      next
    }
    
    # Calculate combinations in chunks to save memory
    chunk_size <- min(10000, choose(length(vars), i))
    total_combs <- choose(length(vars), i)
    
    if(verbose) {
      cat("  Total combinations:", total_combs, "\n")
      cat("  Using chunk size:", chunk_size, "\n")
    }
    
    # Track progress for this feature size
    models_saved_for_size <- 0
    
    # Process combinations in chunks
    for(chunk_start in seq(1, total_combs, by = chunk_size)) {
      chunk_end <- min(chunk_start + chunk_size - 1, total_combs)
      
      # Generate combinations for this chunk
      current_chunk <- utils::combn(vars, i, simplify = FALSE)[chunk_start:chunk_end]
      
      # Track formulas that pass correlation filter
      valid_formulas <- character()
      
      # Apply correlation filtering to each combination
      for(vars_combo in current_chunk) {
        if(length(vars_combo) > 1) {
          # Check correlation between variables
          clean_vars <- gsub("`", "", vars_combo)
          var_indices <- match(clean_vars, colnames(cor_matrix))
          var_cor <- cor_matrix[var_indices, var_indices]
          
          # Skip if any pair has correlation above threshold
          if(any(abs(var_cor[upper.tri(var_cor)]) > cor.threshold)) {
            next
          }
        }
        
        # Add to valid formulas
        valid_formulas <- c(valid_formulas, paste(output, "~", paste(vars_combo, collapse = " + ")))
      }
      
      if(length(valid_formulas) == 0) {
        next
      }
      
      # Calculate R² in parallel for valid formulas
      r_squared_results <- parallel::mclapply(
        valid_formulas,
        function(formula) {
          result <- try(summary(lm(formula, data = data))$r.squared, silent = TRUE)
          if(inherits(result, "try-error")) {
            return(NA)
          }
          return(result)
        },
        mc.cores = max(1, parallel::detectCores() - 1)
      )
      
      # Filter out failed models and create results table
      valid_indices <- !sapply(r_squared_results, is.na)
      valid_formulas <- valid_formulas[valid_indices]
      valid_r2 <- unlist(r_squared_results[valid_indices])
      
      if(length(valid_r2) == 0) {
        next
      }
      
      total_models_evaluated <- total_models_evaluated + length(valid_r2)
      
      # Create results dataframe
      results_df <- data.frame(
        formula = valid_formulas,
        R.sq = valid_r2,
        stringsAsFactors = FALSE
      )
      
      # Save all model results to temporary file
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S_%OS3")
      file_name <- file.path(results_dir, paste0(basename(results_name), "_", i, "_", timestamp, ".csv"))
      write.csv(results_df, file = file_name, row.names = FALSE)
      models_saved_for_size <- models_saved_for_size + nrow(results_df)
      models_saved <- models_saved + nrow(results_df)
    }
    
    if(verbose) {
      cat("  Saved", models_saved_for_size, "models for feature size", i, "\n")
    }
  }
  
  if(verbose) {
    cat("Total models evaluated:", total_models_evaluated, "\n")
    cat("Models saved to temporary files:", models_saved, "\n")
  }
  
  # Get list of all result files from the correct directory
  file_list <- list.files(
    path = results_dir, 
    pattern = paste0('^', basename(results_name), ".*\\.csv$"), 
    full.names = TRUE
  )
  
  if(length(file_list) == 0) {
    warning("No model result files found in directory:", results_dir)
    return(data.frame(
      formula = character(),
      R.sq = numeric(),
      Q.sq = numeric(),
      MAE = numeric(),
      Model = integer(),
      stringsAsFactors = FALSE
    ))
  }
  
  if(verbose) {
    cat("Found", length(file_list), "temporary files with model results\n")
    cat("Finding models with R² ≥", cutoff, "\n")
  }
  
  # First attempt: read files and collect models above initial cutoff
  all_models <- data.frame(formula = character(), R.sq = numeric(), stringsAsFactors = FALSE)
  r_squared_values <- numeric()
  
  for(file in file_list) {
    tryCatch({
      # Read the file
      temp_df <- utils::read.csv(file, stringsAsFactors = FALSE)
      
      # Store all R² values for potential cutoff adjustment
      if(ncol(temp_df) >= 2) {
        r_squared_values <- c(r_squared_values, temp_df[[2]])
        
        # Filter models above cutoff
        models_above_cutoff <- temp_df[temp_df[[2]] >= cutoff, , drop = FALSE]
        if(nrow(models_above_cutoff) > 0) {
          all_models <- rbind(all_models, models_above_cutoff)
        }
      }
    }, error = function(e) {
      warning("Error reading file: ", file, " - ", e$message)
    })
  }
  
  # If not enough models, adaptively adjust cutoff
  if(nrow(all_models) < min_models && length(r_squared_values) > 0) {
    adapted_cutoff <- cutoff
    
    while(nrow(all_models) < min_models && adapted_cutoff > 0.1) {
      # Lower the cutoff
      adapted_cutoff <- adapted_cutoff - 0.05
      
      if(verbose) {
        cat("Adjusting cutoff to", adapted_cutoff, "to find more models\n")
      }
      
      # Clear previous models
      all_models <- data.frame(formula = character(), R.sq = numeric(), stringsAsFactors = FALSE)
      
      # Re-read files with new cutoff
      for(file in file_list) {
        tryCatch({
          temp_df <- utils::read.csv(file, stringsAsFactors = FALSE)
          if(ncol(temp_df) >= 2) {
            models_above_cutoff <- temp_df[temp_df[[2]] >= adapted_cutoff, , drop = FALSE]
            if(nrow(models_above_cutoff) > 0) {
              all_models <- rbind(all_models, models_above_cutoff)
            }
          }
        }, error = function(e) {
          warning("Error reading file: ", file, " - ", e$message)
        })
      }
    }
    
    # If still not enough models, just take the top ones
    if(nrow(all_models) < min_models && length(r_squared_values) >= min_models) {
      if(verbose) {
        cat("Taking top", min_models, "models by R²\n")
      }
      
      cutoff_value <- sort(r_squared_values, decreasing = TRUE)[min(min_models, length(r_squared_values))]
      
      # Clear previous models
      all_models <- data.frame(formula = character(), R.sq = numeric(), stringsAsFactors = FALSE)
      
      # Re-read files with cutoff based on top models
      for(file in file_list) {
        tryCatch({
          temp_df <- utils::read.csv(file, stringsAsFactors = FALSE)
          if(ncol(temp_df) >= 2) {
            models_above_cutoff <- temp_df[temp_df[[2]] >= cutoff_value, , drop = FALSE]
            if(nrow(models_above_cutoff) > 0) {
              all_models <- rbind(all_models, models_above_cutoff)
            }
          }
        }, error = function(e) {
          warning("Error reading file: ", file, " - ", e$message)
        })
      }
    }
  }
  
  # If still no models, take any available
  if(nrow(all_models) == 0 && length(r_squared_values) > 0) {
    if(verbose) {
      cat("No models above cutoff. Taking all available models.\n")
    }
    
    # Read all models without filtering
    for(file in file_list) {
      tryCatch({
        temp_df <- utils::read.csv(file, stringsAsFactors = FALSE)
        if(ncol(temp_df) >= 2) {
          all_models <- rbind(all_models, temp_df)
        }
      }, error = function(e) {
        warning("Error reading file: ", file, " - ", e$message)
      })
    }
  }
  
  # If still no models, handle the error
  if(nrow(all_models) == 0) {
    warning("No valid models found despite adaptive cutoff.")
    
    # Clean up temporary files
    if(verbose) cat("Cleaning up temporary files...\n")
    for(file in file_list) unlink(file)
    
    return(data.frame(
      formula = character(),
      R.sq = numeric(),
      Q.sq = numeric(),
      MAE = numeric(),
      Model = integer(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort models by R²
  all_models <- all_models[order(all_models[[2]], decreasing = TRUE), , drop = FALSE]
  names(all_models) <- c("formula", "R.sq")
  
  # Keep only top models for cross-validation
  if(nrow(all_models) > min_models) {
    all_models <- all_models[1:min_models, , drop = FALSE]
  }
  
  if(verbose) {
    cat("Running cross-validation on", nrow(all_models), "models...\n")
  }
  
  # Perform cross-validation on selected models
  q2.list <- list()
  mae.list <- list()
  
  for(i in 1:nrow(all_models)) {
    if(verbose && i %% 10 == 0) {
      cat("  Cross-validating model", i, "of", nrow(all_models), "\n")
    }
    
    stts <- model.cv(
      formula = all_models$formula[i],
      data = data,
      out.col = out.col,
      folds = folds,
      iterations = iterations
    )
    q2.list[[i]] <- stts[[2]]
    mae.list[[i]] <- stts[[1]]
  }
  
  # Add Q-squared and MAE to the final results
  all_models$Q.sq <- unlist(q2.list)
  all_models$MAE <- unlist(mae.list)
  
  # Sort by Q-squared and add model numbers
  final_models <- all_models[order(all_models$Q.sq, decreasing = TRUE), , drop = FALSE]
  final_models$Model <- seq_len(nrow(final_models))
  
  # Clean up temporary files
  if(verbose) cat("Cleaning up temporary files...\n")
  for(file in file_list) unlink(file)
  
  # Return top 15 models or all models if less than 15
  if(nrow(final_models) > 15) {
    return(final_models[1:15, ])
  } else {
    return(final_models)
  }
}

#' Parallelized Subset Feature Selection with Correlation Threshold (Windows-Compatible)
#'
#' This function provides a Windows-compatible implementation of model.subset.parallel. 
#' It generates all possible subsets of features in a user-defined range (min to max), 
#' ensuring that no two variables in a subset have a Pearson correlation coefficient above the
#' specified correlation threshold. It performs parallel computation using Windows-friendly 
#' approaches (parLapply instead of mclapply) and handles large feature spaces efficiently.
#'
#' @param data A data frame containing the features and the target column.
#' @param out.col An integer indicating the column index of the target variable (default: last column).
#' @param min The minimum number of features to include in a subset (default: 2).
#' @param max The maximum number of features to include in a subset (default: floor(nrow(data) / 5)).
#' @param results_name Base name for temporary files storing model results (default: 'results_df').
#' @param folds Number of folds for cross-validation (default: number of rows in `data`).
#' @param iterations Number of iterations for cross-validation (default: 1).
#' @param cutoff Initial R-squared cutoff for model selection (default: 0.65). The function adaptively adjusts
#'        this value if too few models meet the initial threshold.
#' @param cor.threshold Pearson correlation threshold for feature selection (default: 0.7). Variable pairs with
#'        absolute correlation exceeding this value will not appear together in the same model.
#' @param min_models Minimum number of models to retain for cross-validation (default: 50). The function
#'        adaptively adjusts the cutoff to ensure at least this many models are evaluated.
#' @param verbose Logical indicating whether to print progress information (default: TRUE).
#'
#' @return A data frame with the best models (up to 15), ranked by cross-validated Q-squared value. 
#'         The data frame contains columns for each model's formula, R-squared, Q-squared, MAE, and 
#'         a Model number.
#'
#' @details
#' This Windows-compatible version works identically to model.subset.parallel, but it uses 
#' parallel::makeCluster and parallel::parLapply instead of parallel::mclapply for parallel 
#' processing. This avoids the limitation of mclapply on Windows, which doesn't support 
#' true parallelism and falls back to sequential processing.
#'
#' The function works in several stages:
#' 1. For each feature subset size from `min` to `max`, it generates combinations and filters them
#'    based on the correlation threshold.
#' 2. It evaluates models for these combinations in parallel using Windows-friendly parallelization,
#'    saving results to temporary files.
#' 3. If too few models meet the initial R-squared cutoff, it adaptively lowers the threshold.
#' 4. It performs cross-validation on the top models to compute Q-squared and MAE values.
#' 5. It ranks the final models by Q-squared and returns the top 15.
#'
#' Like its non-Windows counterpart, this function handles memory constraints by processing 
#' combinations in chunks and writing intermediate results to disk.
#'
#' @importFrom stringr str_c
#' @importFrom parallel makeCluster stopCluster clusterExport parLapply detectCores
#' @importFrom utils write.csv read.csv combn
#' @importFrom stats cor lm
#'
#' @export

model.subset.parallel.windows <- function(data, out.col = dim(data)[2],
                                          min = 2, max = floor(dim(data)[1] / 5),
                                          results_name = 'results_df',
                                          folds = nrow(data), iterations = 1,
                                          cutoff = 0.65, cor.threshold = 1,
                                          min_models = 50, verbose = TRUE) {
  
  # Get the output variable and feature names
  output <- stringr::str_c("`", names(data[out.col]), "`")
  vars <- names(data[, -out.col])
  for (i in 1:length(vars)) {
    vars[i] <- stringr::str_c("`", vars[i], "`")
  }
  
  # Compute the correlation matrix for features
  cor_matrix <- cor(data[, -out.col])
  
  # Extract directory path from results_name
  results_dir <- dirname(results_name)
  if(results_dir == ".") results_dir <- getwd()
  
  # Clean up any existing result files with the same name pattern in this directory
  file_pattern <- paste0('^', basename(results_name), ".*\\.csv$")
  old_files <- list.files(path = results_dir, pattern = file_pattern, full.names = TRUE)
  if(length(old_files) > 0) {
    if(verbose) cat("Cleaning up", length(old_files), "old result files...\n")
    for(file in old_files) {
      unlink(file)
    }
  }
  
  # Track total models evaluated
  total_models_evaluated <- 0
  models_saved <- 0
  
  # Process each feature size in range
  for (i in min:max) {
    if(verbose) cat("Processing feature size:", i, "\n")
    
    # Safety check for combinations
    if(choose(length(vars), i) <= 0 || i > length(vars)) {
      if(verbose) cat("  Skipping - invalid feature size\n")
      next
    }
    
    # Calculate combinations in chunks to save memory
    chunk_size <- min(10000, choose(length(vars), i))
    total_combs <- choose(length(vars), i)
    
    if(verbose) {
      cat("  Total combinations:", total_combs, "\n")
      cat("  Using chunk size:", chunk_size, "\n")
    }
    
    # Track progress for this feature size
    models_saved_for_size <- 0
    
    # Process combinations in chunks
    for(chunk_start in seq(1, total_combs, by = chunk_size)) {
      chunk_end <- min(chunk_start + chunk_size - 1, total_combs)
      
      # Generate combinations for this chunk
      current_chunk <- utils::combn(vars, i, simplify = FALSE)[chunk_start:chunk_end]
      
      # Track formulas that pass correlation filter
      valid_formulas <- character()
      
      # Apply correlation filtering to each combination
      for(vars_combo in current_chunk) {
        if(length(vars_combo) > 1) {
          # Check correlation between variables
          clean_vars <- gsub("`", "", vars_combo)
          var_indices <- match(clean_vars, colnames(cor_matrix))
          var_cor <- cor_matrix[var_indices, var_indices]
          
          # Skip if any pair has correlation above threshold
          if(any(abs(var_cor[upper.tri(var_cor)]) > cor.threshold)) {
            next
          }
        }
        
        # Add to valid formulas
        valid_formulas <- c(valid_formulas, paste(output, "~", paste(vars_combo, collapse = " + ")))
      }
      
      if(length(valid_formulas) == 0) {
        next
      }
      
      # WINDOWS MODIFICATION: Use parLapply instead of mclapply
      # Create a cluster for parallel processing
      n_cores <- max(1, parallel::detectCores() - 1)
      cl <- parallel::makeCluster(n_cores)
      
      # Export necessary objects to the cluster
      parallel::clusterExport(cl, c("data"), envir = environment())
      
      # Calculate R² in parallel for valid formulas
      r_squared_results <- parallel::parLapply(
        cl,
        valid_formulas,
        function(formula) {
          result <- try(summary(lm(formula, data = data))$r.squared, silent = TRUE)
          if(inherits(result, "try-error")) {
            return(NA)
          }
          return(result)
        }
      )
      
      # Close the cluster
      parallel::stopCluster(cl)
      
      # Filter out failed models and create results table
      valid_indices <- !sapply(r_squared_results, is.na)
      valid_formulas <- valid_formulas[valid_indices]
      valid_r2 <- unlist(r_squared_results[valid_indices])
      
      if(length(valid_r2) == 0) {
        next
      }
      
      total_models_evaluated <- total_models_evaluated + length(valid_r2)
      
      # Create results dataframe
      results_df <- data.frame(
        formula = valid_formulas,
        R.sq = valid_r2,
        stringsAsFactors = FALSE
      )
      
      # Save all model results to temporary file
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S_%OS3")
      file_name <- file.path(results_dir, paste0(basename(results_name), "_", i, "_", timestamp, ".csv"))
      write.csv(results_df, file = file_name, row.names = FALSE)
      models_saved_for_size <- models_saved_for_size + nrow(results_df)
      models_saved <- models_saved + nrow(results_df)
    }
    
    if(verbose) {
      cat("  Saved", models_saved_for_size, "models for feature size", i, "\n")
    }
  }
  
  if(verbose) {
    cat("Total models evaluated:", total_models_evaluated, "\n")
    cat("Models saved to temporary files:", models_saved, "\n")
  }
  
  # Get list of all result files from the correct directory
  file_list <- list.files(
    path = results_dir, 
    pattern = paste0('^', basename(results_name), ".*\\.csv$"), 
    full.names = TRUE
  )
  
  if(length(file_list) == 0) {
    warning("No model result files found in directory:", results_dir)
    return(data.frame(
      formula = character(),
      R.sq = numeric(),
      Q.sq = numeric(),
      MAE = numeric(),
      Model = integer(),
      stringsAsFactors = FALSE
    ))
  }
  
  if(verbose) {
    cat("Found", length(file_list), "temporary files with model results\n")
    cat("Finding models with R² ≥", cutoff, "\n")
  }
  
  # First attempt: read files and collect models above initial cutoff
  all_models <- data.frame(formula = character(), R.sq = numeric(), stringsAsFactors = FALSE)
  r_squared_values <- numeric()
  
  for(file in file_list) {
    tryCatch({
      # Read the file
      temp_df <- utils::read.csv(file, stringsAsFactors = FALSE)
      
      # Store all R² values for potential cutoff adjustment
      if(ncol(temp_df) >= 2) {
        r_squared_values <- c(r_squared_values, temp_df[[2]])
        
        # Filter models above cutoff
        models_above_cutoff <- temp_df[temp_df[[2]] >= cutoff, , drop = FALSE]
        if(nrow(models_above_cutoff) > 0) {
          all_models <- rbind(all_models, models_above_cutoff)
        }
      }
    }, error = function(e) {
      warning("Error reading file: ", file, " - ", e$message)
    })
  }
  
  # If not enough models, adaptively adjust cutoff
  if(nrow(all_models) < min_models && length(r_squared_values) > 0) {
    adapted_cutoff <- cutoff
    
    while(nrow(all_models) < min_models && adapted_cutoff > 0.1) {
      # Lower the cutoff
      adapted_cutoff <- adapted_cutoff - 0.05
      
      if(verbose) {
        cat("Adjusting cutoff to", adapted_cutoff, "to find more models\n")
      }
      
      # Clear previous models
      all_models <- data.frame(formula = character(), R.sq = numeric(), stringsAsFactors = FALSE)
      
      # Re-read files with new cutoff
      for(file in file_list) {
        tryCatch({
          temp_df <- utils::read.csv(file, stringsAsFactors = FALSE)
          if(ncol(temp_df) >= 2) {
            models_above_cutoff <- temp_df[temp_df[[2]] >= adapted_cutoff, , drop = FALSE]
            if(nrow(models_above_cutoff) > 0) {
              all_models <- rbind(all_models, models_above_cutoff)
            }
          }
        }, error = function(e) {
          warning("Error reading file: ", file, " - ", e$message)
        })
      }
    }
    
    # If still not enough models, just take the top ones
    if(nrow(all_models) < min_models && length(r_squared_values) >= min_models) {
      if(verbose) {
        cat("Taking top", min_models, "models by R²\n")
      }
      
      cutoff_value <- sort(r_squared_values, decreasing = TRUE)[min(min_models, length(r_squared_values))]
      
      # Clear previous models
      all_models <- data.frame(formula = character(), R.sq = numeric(), stringsAsFactors = FALSE)
      
      # Re-read files with cutoff based on top models
      for(file in file_list) {
        tryCatch({
          temp_df <- utils::read.csv(file, stringsAsFactors = FALSE)
          if(ncol(temp_df) >= 2) {
            models_above_cutoff <- temp_df[temp_df[[2]] >= cutoff_value, , drop = FALSE]
            if(nrow(models_above_cutoff) > 0) {
              all_models <- rbind(all_models, models_above_cutoff)
            }
          }
        }, error = function(e) {
          warning("Error reading file: ", file, " - ", e$message)
        })
      }
    }
  }
  
  # If still no models, take any available
  if(nrow(all_models) == 0 && length(r_squared_values) > 0) {
    if(verbose) {
      cat("No models above cutoff. Taking all available models.\n")
    }
    
    # Read all models without filtering
    for(file in file_list) {
      tryCatch({
        temp_df <- utils::read.csv(file, stringsAsFactors = FALSE)
        if(ncol(temp_df) >= 2) {
          all_models <- rbind(all_models, temp_df)
        }
      }, error = function(e) {
        warning("Error reading file: ", file, " - ", e$message)
      })
    }
  }
  
  # If still no models, handle the error
  if(nrow(all_models) == 0) {
    warning("No valid models found despite adaptive cutoff.")
    
    # Clean up temporary files
    if(verbose) cat("Cleaning up temporary files...\n")
    for(file in file_list) unlink(file)
    
    return(data.frame(
      formula = character(),
      R.sq = numeric(),
      Q.sq = numeric(),
      MAE = numeric(),
      Model = integer(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort models by R²
  all_models <- all_models[order(all_models[[2]], decreasing = TRUE), , drop = FALSE]
  names(all_models) <- c("formula", "R.sq")
  
  # Keep only top models for cross-validation
  if(nrow(all_models) > min_models) {
    all_models <- all_models[1:min_models, , drop = FALSE]
  }
  
  if(verbose) {
    cat("Running cross-validation on", nrow(all_models), "models...\n")
  }
  
  # Perform cross-validation on selected models
  q2.list <- list()
  mae.list <- list()
  
  for(i in 1:nrow(all_models)) {
    if(verbose && i %% 10 == 0) {
      cat("  Cross-validating model", i, "of", nrow(all_models), "\n")
    }
    
    # WINDOWS MODIFICATION: Use Windows-specific version of model.cv.parallel
    stts <- model.cv(
      formula = all_models$formula[i],
      data = data,
      out.col = out.col,
      folds = folds,
      iterations = iterations
    )
    q2.list[[i]] <- stts[[2]]
    mae.list[[i]] <- stts[[1]]
  }
  
  # Add Q-squared and MAE to the final results
  all_models$Q.sq <- unlist(q2.list)
  all_models$MAE <- unlist(mae.list)
  
  # Sort by Q-squared and add model numbers
  final_models <- all_models[order(all_models$Q.sq, decreasing = TRUE), , drop = FALSE]
  final_models$Model <- seq_len(nrow(final_models))
  
  # Clean up temporary files
  if(verbose) cat("Cleaning up temporary files...\n")
  for(file in file_list) unlink(file)
  
  # Return top 15 models or all models if less than 15
  if(nrow(final_models) > 15) {
    return(final_models[1:15, ])
  } else {
    return(final_models)
  }
}

####### ----------------------------------------------------#####
####### -------------------User Functions-------------------#####
####### ----------------------------------------------------#####

#' Generate Linear Models from Feature Subset Selection with Parallel Processing
#'
#' This function performs feature subset selection to identify optimal linear models and saves
#' the results to a timestamped directory. It handles data preprocessing, cross-validation,
#' and optional prediction on left-out samples.
#'
#' @param dataset Path to a CSV file containing the data. First column should contain row identifiers,
#'        and the last column should contain the outcome variable.
#' @param min Minimum number of features to include in a subset (default: 2).
#' @param max Maximum number of features to include in a subset (default: NULL, which automatically
#'        sets the value to min(floor(nrow(data)/5), 5)).
#' @param leave.out Vector of row identifiers to exclude from model training. These samples
#'        will be used for external validation if specified (default: '').
#' @param folds Number of folds for cross-validation (default: number of rows in the dataset).
#' @param iterations Number of iterations for cross-validation (default: 1).
#' @param cutoff Initial R-squared cutoff for model selection (default: 0.50). The function adaptively
#'        adjusts this value if too few models meet the initial threshold.
#' @param cor.threshold Pearson correlation threshold for feature selection (default: 1.0). Variable pairs
#'        with absolute correlation exceeding this value will not appear together in the same model.
#' @param min_models Minimum number of models to retain for cross-validation (default: 50). The function
#'        adaptively adjusts the cutoff to ensure at least this many models are evaluated.
#' @param verbose Logical indicating whether to print progress information (default: TRUE).
#'
#' @return A data frame with the best models (up to 15), ranked by cross-validated Q-squared value.
#'         The data frame contains columns for each model's formula, R-squared, Q-squared, MAE, and
#'         a Model number. The results are also saved to a CSV file in a timestamped directory.
#'
#' @details
#' The function works in several stages:
#' 1. It creates a timestamped directory for storing results
#' 2. It reads and preprocesses the dataset, including scaling predictors
#' 3. It separates any leave-out samples for later prediction
#' 4. It runs model.subset.parallel to identify optimal models
#' 5. It saves the results to a CSV file in the timestamped directory
#' 6. If leave-out samples were provided, it makes predictions using the best model
#'    and saves these to a separate CSV file
#'
#' The function handles error conditions gracefully and provides detailed progress information
#' when verbose=TRUE. It automatically determines appropriate parameter values when not explicitly
#' provided, making it suitable for both exploratory analysis and automated workflows.
#'
#' @importFrom utils read.csv write.csv
#' @importFrom stats complete.cases as.formula predict
#' @importFrom tools file_path_sans_ext
#' @importFrom knitr kable
#' 
#' @aliases models.list.parallel
#' @export
models.list.parallel <- function(dataset,
                                 min = 2,
                                 max = NULL,
                                 leave.out = '',
                                 folds = nrow(read.csv(dataset)), 
                                 iterations = 1,
                                 cutoff = 0.50,
                                 cor.threshold = 1.0, # Allow any correlation by default
                                 min_models = 50,
                                 verbose = TRUE) {
  
  # Create directory name based on parameters and timestamp
  dir_name <- paste0(
    tools::file_path_sans_ext(basename(dataset)),
    '_min', as.character(min),
    '_max', ifelse(is.null(max), "auto", as.character(max)),
    '_', format(Sys.time(), "%Y%m%d_%H%M%S")
  )
  
  # Create directory if it doesn't exist
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  # Get absolute path to ensure consistent file handling
  dir_path <- normalizePath(dir_name)
  
  # Modify results_name to include directory path
  results_file_path <- file.path(dir_path, "results_df")
  
  if(verbose) {
    cat("Reading and processing dataset:", dataset, "\n")
    cat("Results will be saved in directory:", dir_name, "\n")
  }
  
  # Read and preprocess data
  mod_data <- tryCatch({
    utils::read.csv(dataset, stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    stop("Error reading dataset: ", e$message)
  })
  
  if(ncol(mod_data) <= 1) {
    stop("Dataset must have at least two columns (row names and one variable)")
  }
  
  # Store row names (first column) and remove it
  RN <- mod_data[,1]
  mod_data <- mod_data[,-1]
  
  # Remove rows with missing values
  mod_data <- mod_data[complete.cases(mod_data), ]
  
  if(nrow(mod_data) == 0) {
    stop("No complete cases found in the dataset after removing missing values")
  }
  
  CN <- names(mod_data)
  
  # Scale predictors (all columns except the last one)
  predictors <- mod_data[, 1:(ncol(mod_data) - 1), drop = FALSE]
  outcome <- mod_data[, ncol(mod_data)]
  
  # Scale remaining predictors
  if(ncol(predictors) > 0) {
    scaled_predictors <- scale(predictors, center = TRUE, scale = TRUE)
    
    # Check for NaN or Inf values after scaling
    if(any(is.na(scaled_predictors)) || any(is.infinite(scaled_predictors))) {
      warning("Scaling produced NaN or Inf values. Using original unscaled predictors.")
      scaled_predictors <- as.data.frame(predictors)
    } else {
      scaled_predictors <- as.data.frame(scaled_predictors)
    }
    
    # Combine scaled predictors with outcome
    mod_data <- data.frame(scaled_predictors, outcome = outcome)
    names(mod_data) <- c(names(scaled_predictors), CN[length(CN)])
  } else {
    warning("No valid predictor columns remain.")
    return(data.frame())
  }
  
  row.names(mod_data) <- RN
  
  # Extract observations to leave out
  if(length(leave.out) > 0 && any(leave.out != '')) {
    leave_out_indices <- row.names(mod_data) %in% leave.out
    if(any(leave_out_indices)) {
      pred.data <- mod_data[leave_out_indices, , drop = FALSE]
      mod_data <- mod_data[!leave_out_indices, , drop = FALSE]
      
      if(verbose) {
        cat("Left out", nrow(pred.data), "samples for testing\n")
      }
    } else {
      warning("Specified leave.out samples not found in dataset")
      pred.data <- NULL
    }
  } else {
    pred.data <- NULL
  }
  
  # Set default max if not provided
  if (is.null(max)) {
    max <- min(floor(nrow(mod_data) / 5), 5)  # Cap at 5 for small datasets
    if(verbose) cat("Setting max number of features to:", max, "\n")
  }
  
  # Default number of folds - cap at number of observations
  folds <- min(folds, nrow(mod_data))
  if(verbose) cat("Setting number of cross-validation folds to:", folds, "\n")
  
  if(verbose) {
    cat("Running model subset selection with parameters:\n")
    cat("- Min features:", min, "\n")
    cat("- Max features:", max, "\n")
    cat("- Initial R² cutoff:", cutoff, "\n")
    cat("- Correlation threshold:", cor.threshold, "\n")
    cat("- Min models for CV:", min_models, "\n")
    cat("- CV folds:", folds, "\n")
    cat("- CV iterations:", iterations, "\n")
    cat("- Number of observations:", nrow(mod_data), "\n")
    cat("- Number of variables:", ncol(mod_data) - 1, "\n")
  }
  
  # Check if there are enough observations and variables
  if(nrow(mod_data) < folds) {
    warning("Number of observations (", nrow(mod_data), ") is less than number of folds (", folds, ").",
            "Setting folds to ", nrow(mod_data))
    folds <- nrow(mod_data)
  }
  
  if(ncol(mod_data) - 1 < min) {
    stop("Too few variables (", ncol(mod_data) - 1, ") to meet minimum feature requirement (", min, ")")
  }
  
  # Run model.subset.parallel with robust error handling
  models <- tryCatch({
    model.subset.parallel(
      data = mod_data,
      out.col = ncol(mod_data),  # Last column
      min = min,
      max = max,
      results_name = results_file_path,
      folds = folds, 
      iterations = iterations,
      cutoff = cutoff,
      cor.threshold = cor.threshold,
      min_models = min_models,
      verbose = verbose
    )
  }, error = function(e) {
    cat("Error in model.subset.parallel:", e$message, "\n")
    
    # Return empty data frame with appropriate structure
    data.frame(
      formula = character(),
      R.sq = numeric(),
      Q.sq = numeric(),
      MAE = numeric(),
      Model = integer(),
      stringsAsFactors = FALSE
    )
  })
  
  if(nrow(models) == 0) {
    warning("No models found. Check your data and parameters.")
  } else {
    if(verbose) {
      cat("\nFound", nrow(models), "models.\n")
      cat("Top", min(5, nrow(models)), "models by Q²:\n")
      print(knitr::kable(head(models, 5)))
    }
    
    # Save models list in the results directory
    output_file <- file.path(
      dir_path,
      paste0(
        tools::file_path_sans_ext(basename(dataset)),
        '_models_list.csv'
      )
    )
    
    utils::write.csv(models, output_file, row.names = FALSE)
    
    if(verbose) {
      cat("\nResults saved to:", output_file, "\n")
    }
    
    # If there are left-out samples and at least one model, make predictions
    if(!is.null(pred.data) && nrow(pred.data) > 0 && nrow(models) > 0) {
      best_model <- tryCatch({
        lm(as.formula(models$formula[1]), data = mod_data)
      }, error = function(e) {
        warning("Error fitting best model for predictions: ", e$message)
        NULL
      })
      
      if(!is.null(best_model)) {
        predictions <- tryCatch({
          predict(best_model, newdata = pred.data)
        }, error = function(e) {
          warning("Error making predictions: ", e$message)
          NULL
        })
        
        if(!is.null(predictions)) {
          actual_values <- pred.data[[ncol(pred.data)]]
          
          prediction_results <- data.frame(
            Sample = row.names(pred.data),
            Actual = actual_values,
            Predicted = predictions,
            Error = actual_values - predictions
          )
          
          # Save predictions
          pred_file <- file.path(
            dir_path,
            paste0(
              tools::file_path_sans_ext(basename(dataset)),
              '_predictions.csv'
            )
          )
          
          utils::write.csv(prediction_results, pred_file, row.names = FALSE)
          
          if(verbose) {
            cat("\nPredictions for left-out samples saved to:", pred_file, "\n")
          }
        }
      }
    }
  }
  
  return(models)
}

#'  Generate models, cross validate, and create visualizations with parallel processing
#'
#' Screen for models using parallel processing, cross validate them and plot. Designed for interactive work.
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
model.report.parallel <- function(dataset, min = 2, max = NULL,
                                  leave.out = '', predict = FALSE,
                                  what.model = NULL, cor.threshold = 0.7,
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
  
  # Run model search with correlation threshold - USING PARALLEL VERSION
  cat("\nRunning parallel model search with correlation threshold:", cor.threshold, "\n")
  tryCatch({
    # The only real change is here - using model.subset.parallel instead of model.subset
    models <- model.subset.parallel(
      data = mod_data, 
      out.col = ncol(mod_data),
      min = min, 
      max = max, 
      iterations = 1, 
      cutoff = cutoff, 
      cor.threshold = cor.threshold
    )
    
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


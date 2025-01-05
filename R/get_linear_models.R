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
                         cutoff = 0.85, cor.threshold = 0.7) {
  
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
                           newdata = data,
                           interval = 'pre',
                           level = 0.9)
  plot.dat <- data.frame(cbind(data$output, pred_interval))
  colnames(plot.dat) <- c('Measured', 'Predicted', 'lwr', 'upr')
  rownames(plot.dat) <- row.names(data)

  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o_",'2-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"m_",'3-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"p_",'4-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o4-",'2,4-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"m3-",'3,5-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o3-",'2,3-')
  

  plot.dat <- dplyr::mutate(plot.dat, Position = rep(NA, nrow(plot.dat)))
  for (i in 1:nrow(plot.dat)) {
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

  plot <- ggplot2::ggplot(plot.dat, ggplot2::aes(x = Measured, y = Predicted)) +
    ggplot2::geom_point(size = 2, shape = 15, ggplot2::aes(color = Position)) +
    ggplot2::stat_smooth(ggplot2::aes(y = lwr), color = "cadetblue", linetype = "dashed",
                se = F, method = 'lm', fullrange = T, linewidth = 0.8) +
    ggplot2::stat_smooth(ggplot2::aes(y = upr), color = "cadetblue", linetype = "dashed",
                se = F, method = 'lm', fullrange = T, linewidth = 0.8) +
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
    ggplot2::scale_color_manual(values = c(Ph = "black", meta = 'tan1',
                                  para = '#66a182',ortho = '#d1495b')) +
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
    ggplot2::theme(text = ggplot2::element_text(family = 'Arial'))
  plot
}

####### ----------------------------------------------------#####
####### -------------------User Functions-------------------#####
####### ----------------------------------------------------#####

#'  Generate models and a plot
#'
#' Screen for models, cross validate them and plot. Designed for interactive work.
#' @param dataset a dataframe with outcome column (must be named 'output')
#' @param min minimum # of features (default = 2)
#' @param max max # of features (defaults = # of observations / 5)
#' @param leave.out name of observations to leave out (e.g. 'p_Br')
#' @param predict if leave.out is not empty, should a prediction
#' of it be computed?
#' @param what.model if not in an interactive session, what model should be used?
#' @return models list, CV results for k=3/5/LOO and a plot of choice. interactive.
#' @export
model.report <- function(dataset, min = 2, max = floor(dim(mod_data)[1] / 5),
                          leave.out = '', predict = F,
                         what.model = NULL) {
  default::default(data.frame) <- list(check.names = FALSE)
  cat(tools::file_path_sans_ext(basename(dataset)))
  mod_data <- data.frame(data.table::fread(dataset, header = T,
                                           check.names = F))
  RN <- mod_data[,1]
  mod_data <- mod_data[,-1]
  mod_data <- mod_data[complete.cases(mod_data), ]
  CN <- names(mod_data)
  mod_data <- data.frame(cbind(scale(mod_data[,1:dim(mod_data)[2] - 1], T, T), mod_data$output))
  names(mod_data) <- CN
  row.names(mod_data) <- RN
  pred.data <- mod_data[row.names(mod_data) %in% leave.out, ]
  mod_data <- mod_data[!(row.names(mod_data) %in% leave.out), ]
  cutoff <- 0.85
  models <- model.subset(data = mod_data, min = min, max = max, iterations = 1, cutoff = cutoff, leave.out = leave.out)
  tab <- knitr::kable(models)
  print(tab)
  if (is.null(what.model)) what.model <- readline('Choose the model you would like to plot (line number): ')
  if (is.character(what.model)) what.model <- as.numeric(what.model)
  mod.sum <- summary(lm(models$formula[what.model], mod_data))$coefficients
  cat('
  Model Coefficients')
  colnames(mod.sum)[4] <- 'p value'
  k.mod <- knitr::kable(mod.sum)
  print(k.mod)
  cv_3fold <- model.cv(models[what.model,1], mod_data, dim(mod_data)[2], 3, 50)
  dt3 <- data.frame(cv_3fold[[2]], cv_3fold[[1]])
  names(dt3) <- c('Q2', 'MAE')
  cat('
  3-fold CV')
  tab_dt3 <- knitr::kable(dt3)
  print(tab_dt3)
  cv_5fold <- model.cv(models[what.model,1], mod_data, dim(mod_data)[2], 5, 50)
  dt5 <- data.frame(cv_5fold[[2]], cv_5fold[[1]])
  names(dt5) <- c('Q2', 'MAE')
  cat('
  5-fold CV')
  tab_dt5 <- knitr::kable(dt5)
  print(tab_dt5)
  if (predict == T) {
    prediction <- predict(lm(models$formula[what.model], mod_data), pred.data)
    real <- pred.data$output
    prd.tab <- data.frame(prediction, real)
    names(prd.tab) <- c('OOS Pred', 'OOS Measured')
    k.prd.tab <- knitr::kable(prd.tab)
    print(k.prd.tab)
  }
  mod_data_unn <- data.frame(data.table::fread(dataset, header = T))
  mod.sum.unnormalized <- summary(lm(models$formula[what.model], mod_data_unn))$coefficients
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

  model = models$formula[what.model]
  best.mod <- lm(model, data = mod_data)
  pred_interval <- predict(best.mod,
                           newdata = mod_data,
                           interval = 'pre',
                           level = 0.9)
  plot.dat <- data.frame(cbind(mod_data$output, pred_interval))
  colnames(plot.dat) <- c('Measured', 'Predicted', 'lwr', 'upr')
  rownames(plot.dat) <- row.names(mod_data)
  if (predict == T) {
    prd.plot.dat <- data.frame(real, prediction,
                               lwr = rep(NA, nrow(pred.data)),
                               upr = rep(NA, nrow(pred.data))
    )
    colnames(prd.plot.dat) <- c('Measured', 'Predicted', 'lwr', 'upr')
    plot.dat <- data.frame(rbind(plot.dat, prd.plot.dat))
  }
  
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
  if (predict == T) {
    plot.dat <- dplyr::mutate(plot.dat, 
                              shapes = c(rep(18, nrow(mod_data)),
                                         rep(19, nrow(pred.data))))
  } else {
    plot.dat <- dplyr::mutate(plot.dat, 
                              shapes = c(rep(18, nrow(mod_data))))
  }
  plot <- suppressMessages(ggplot2::ggplot(plot.dat, ggplot2::aes(x = Measured, y = Predicted)) +
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
                                            legend.position.inside = c(2,2)) +
                             ggplot2::scale_color_manual('', values = c(Ph = "black", meta = 'tan1',
                                                                        para = '#66a182',ortho = '#d1495b', external = 'steelblue4')) +
                             ggplot2::xlim(min(plot.dat[1:nrow(mod_data),3]), max(plot.dat[1:nrow(mod_data),4])) +
                             ggplot2::ylim(min(plot.dat[1:nrow(mod_data),3]), max(plot.dat[1:nrow(mod_data),4])) +
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
                                               x = min(plot.dat[1:nrow(mod_data),3]),
                                               y = max(plot.dat[1:nrow(mod_data),2]), label = annotations,
                                               parse = F,
                                               hjust = "left", vjust = 0) +
                             ggplot2::ggtitle(equation))
  plot
}

#'  Generate a model report and a plot from a model list (after model.subset)
#'
#' Screen for models, cross validate them and plot. Designed for interactive work.
#' @param dataset a dataframe with outcome column (must be named 'output')
#' @param model.list minimum # of features (default = 2)
#' @param out.col column number of y value
#' @param leave.out name of observations to leave out (e.g. 'p_Br') and retrain the model
#' @param predict if leave.out is not empty, should a prediction
#' of it be computed?
#' @param what.model if not in an interactive session, what model should be used?
#' @return models list, CV results for k=3/5/LOO and a plot of choice. interactive.
#' @export
model.report.from.list <- function(dataset, model.list, out.col = 'output',
                                   leave.out = '', predict = F,
                                   what.model = 1) {
  default::default(data.frame) <- list(check.names = FALSE)
  cat(tools::file_path_sans_ext(basename(dataset)))
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
  pred.data <- mod_data[row.names(mod_data) %in% leave.out, ]
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
  if (predict == T) {
    prediction <- predict(lm(models$formula[what.model], mod_data), pred.data)
    real <- pred.data[, dim(mod_data)[2]]
    prd.tab <- data.frame(prediction, real)
    names(prd.tab) <- c('OOS Pred', 'OOS Measured')
    k.prd.tab <- knitr::kable(prd.tab)
    print(k.prd.tab)
  }
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
  info.table[2,1] <- as.character(round(models$Q.sq[what.model], 2))
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
  
  model = models$formula[what.model]
  best.mod <- lm(model, data = mod_data)
  pred_interval <- predict(best.mod,
                           newdata = mod_data,
                           interval = 'pre',
                           level = 0.9)
  plot.dat <- data.frame(cbind(mod_data[dim(mod_data)[2]], pred_interval))
  colnames(plot.dat) <- c('Measured', 'Predicted', 'lwr', 'upr')
  rownames(plot.dat) <- row.names(mod_data)
  if (predict == T) {
    prd.plot.dat <- data.frame(real, prediction,
                               lwr = rep(NA, nrow(pred.data)),
                               upr = rep(NA, nrow(pred.data))
    )
    colnames(prd.plot.dat) <- c('Measured', 'Predicted', 'lwr', 'upr')
    plot.dat <- data.frame(rbind(plot.dat, prd.plot.dat))
  }
  
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
  if (predict == T) {
    plot.dat <- dplyr::mutate(plot.dat, 
                              shapes = c(rep(18, nrow(mod_data)),
                                         rep(19, nrow(pred.data))))
  } else {
    plot.dat <- dplyr::mutate(plot.dat, 
                              shapes = c(rep(18, nrow(mod_data))))
  }
  plot <- suppressMessages(ggplot2::ggplot(plot.dat, ggplot2::aes(x = Measured, y = Predicted)) +
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
                                            legend.position.inside = c(2,2)) +
                             ggplot2::scale_color_manual('', values = c(Ph = "black", meta = 'tan1',
                                                                        para = '#66a182',ortho = '#d1495b', external = 'steelblue4')) +
                             ggplot2::xlim(min(plot.dat[1:nrow(mod_data),3]), max(plot.dat[1:nrow(mod_data),4])) +
                             ggplot2::ylim(min(plot.dat[1:nrow(mod_data),3]), max(plot.dat[1:nrow(mod_data),4])) +
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
                                               x = min(plot.dat[1:nrow(mod_data),3]),
                                               y = max(plot.dat[1:nrow(mod_data),2]), label = annotations,
                                               parse = F,
                                               hjust = "left", vjust = 0) +
                             ggplot2::ggtitle(equation))
  plot
}


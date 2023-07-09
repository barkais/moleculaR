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
  MAE.3cv <- mean(sq.dis.cv$abs.pred...data...out.col..)
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
#'
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

#' Brute force model search with a min and max number of
#' features.
#'
#' Ranks models based on K-fold CV Q2
#'
#' @param data data frame with output column
#' @param out.col number of output column
#' @param min minimum # of features (default = 2)
#' @param max max # of features (defaults = # of observations / 5)
#' @param folds  defaults to nrow(data)
#' @param iterations defaults to 1 (LOOCV)
#' @param cutoff search for Q2 above 0.85 (if there isn't will look for lower)
#' @param cross.terms if TRUE includes feature interactions (explosive - try avoiding)
#'
#' @return table with 10 best models (at max)
#' @export
model.subset <- function(data, out.col = dim(data)[2],
                      min = 2, max = floor(dim(data)[1] / 5),
                      folds = nrow(data), iterations = 1,
                      cutoff = 0.85, cross.terms = F) {
  output <- stringr::str_c("`", names(data[out.col]), "`")
  vars <- names(data[, -out.col])
  for (i in 1:length(vars)) {
    vars[i] <- stringr::str_c("`", vars[i], "`")
  }
  comb.list <- list()
  ols.list <- list()
  q2.list <- list()
  mae.list <- list()
  if (cross.terms == T) {
    max <- max - 1
  }
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
  if (cross.terms == T) {
    cross.list <- list()
    cross.list[[1]] <- data.frame(aperm(combn(vars, 2)), stringsAsFactors = F)
    cross.list[[2]] <- do.call(paste,
                               c(cross.list[[1]][names(cross.list[[1]])],
                                 sep = " : "))
    repli_single <- do.call(rbind,
                            replicate(length(cross.list[[2]]),
                                      forms, simplify = F))
    repli_single <- repli_single[order(repli_single), ]
    repli_cross <- do.call(c, replicate(nrow(forms),
                                        cross.list[[2]],
                                        simplify = F))
    new_comb <- list()
    for (i in 1:length(repli_cross)) {
      new_comb[[i]] <- paste(repli_single[[i]],
                             repli_cross[[i]], sep = " + ")
    }
    orig.forms <- forms
    forms <- data.frame(do.call(rbind,
                                plyr::compact(new_comb)),
                        stringsAsFactors = F)
    names(forms) <- "formula"
    forms <- data.frame(rbind(orig.forms, forms))
    colnames(forms) <- "formula"
  }
  forms$formula <- stringr::str_c(output, " ~ ", forms$formula)
  if (max > 3) { # a work around for high dimension models
    ols.low <- list()
    forms.low <- forms[stringr::str_count(forms$formula, pattern = "\\+") <= 2, ]
    forms.high <- forms[stringr::str_count(forms$formula, pattern = "\\+") >= (min - 1), ]
    ols.low <- unlist(lapply(forms.low, function(x) summary(lm(x, data = data))$r.squared))
    forms.low <- data.frame(forms.low, ols.low)
    forms.low <- forms.low[order(forms.low$ols.low, decreasing = T),]
    forms.low <- forms.low[1:10, ]
    forms.low$forms.low <- stringr::str_remove(forms.low$forms.low, '`output` ~ ')
    broken.forms <- stringr::str_split(forms.low$forms.low, ' \\+ ')
    new.forms <- list()
    for (i in 1:10) {
      new.forms[[i]] <- forms.high[which(apply(sapply(X = broken.forms[[i]],
                                                      FUN = grepl,
                                                      forms.high),
                                               1,
                                               all))]
    }
    forms <- data.frame(unlist(new.forms))
    names(forms) <- 'formula'
  }
  ols.list <- lapply(forms$formula, function(x) summary(lm(x, data = data))$r.squared)
  forms[, 2] <- do.call(rbind, ols.list)
  names(forms)[2] <- "R.sq"
  forms.cut <- forms[forms$R.sq > cutoff, ]
  while (nrow(forms.cut) == 0) {
    cutoff <- cutoff - 0.1
    forms.cut <- forms[forms$R.sq > cutoff, ]
  }
  if (nrow(forms.cut) >= 10) forms.cut <- forms.cut[1:10, ]
  for (i in 1:dim(forms.cut)[1]) {
    stts <- moleculaR:::model.cv(forms.cut[i, 1], data, out.col, folds, iterations)
    q2.list[[i]] <- stts[2]
    mae.list[[i]] <- stts[1]
  }
  forms.cut[, 3] <- data.table::transpose(do.call(rbind, q2.list))
  forms.cut[, 4] <- data.table::transpose(do.call(rbind, mae.list))
  names(forms.cut)[3:4] <- c("Q.sq", "MAE")
  forms.cut <- dplyr::arrange(forms.cut, desc(forms.cut$Q.sq))
  forms.cut <- dplyr::mutate(forms.cut, Model = seq(1, nrow(forms.cut)))
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
#' @param ext.val if leave.out is not empty, should this list be used as a test?
#' @param what.model if not in an interactive session, what model should be used?
#' @return models list, CV results for k=3/5/LOO and a plot of choice. interactive.
#' @export
model.report <- function(dataset, min = 2, max = floor(dim(mod_data)[1] / 5),
                          leave.out = '', predict = F, ext.val = F,
                         what.model = NULL) {
  cat(tools::file_path_sans_ext(basename(dataset)))
  mod_data <- data.frame(data.table::fread(dataset, header = T,
                                           check.names = T))
  RN <- mod_data[,1]
  mod_data <- mod_data[,-1]
  mod_data <- mod_data[complete.cases(mod_data), ]
  CN <- names(mod_data)
  mod_data <- data.frame(cbind(scale(mod_data[,1:dim(mod_data)[2] - 1], T, T), mod_data$output))
  names(mod_data) <- CN
  row.names(mod_data) <- RN
  pred.data <- mod_data[row.names(mod_data) %in% leave.out, ]
  mod_data <- mod_data[!(row.names(mod_data) %in% leave.out), ]
  ifelse(ext.val == T, cutoff = 0.8, cutoff = 0.9)
  models <- model.subset(mod_data, min = min, max = max, iterations = 1, cutoff = 0.9)
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
  cv_3fold <- moleculaR:::model.cv(models[what.model,1], mod_data, dim(mod_data)[2], 3, 50)
  dt3 <- data.frame(cv_3fold[[2]], cv_3fold[[1]])
  names(dt3) <- c('Q2', 'MAE')
  cat('
  3-fold CV')
  tab_dt3 <- knitr::kable(dt3)
  print(tab_dt3)
  cv_5fold <- moleculaR:::model.cv(models[what.model,1], mod_data, dim(mod_data)[2], 5, 50)
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


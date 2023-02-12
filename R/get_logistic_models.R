#' Stratify data sets by a group
#'
#' User defined stratification by a class column. 
#' Designated for k-fold CV. 
#'
#' @param df a dataframe with a categorial column
#' @param group group by which stratified
#' @param size a numeric above or below 1, with above sampling the number
#' and below sampling the ratio.
#' @keywords internal
#' @return a stratum df
stratified <- function(df, group, size) {
  df.interaction <- interaction(df[group], drop = TRUE)
  df.table <- table(df.interaction)
  df.split <- split(df, df.interaction)
  if (length(size) > 1) {
    if (length(size) != length(df.split))
      stop("Number of groups is ", length(df.split),
           " but number of sizes supplied is ", length(size))
    if (is.null(names(size))) {
      n <- setNames(size, names(df.split))
    } else {
      ifelse(all(names(size) %in% names(df.split)),
             n <- size[names(df.split)],
             stop("Named vector supplied with names ",
                  paste(names(size), collapse = ", "),
                  "\n but the names for the group levels are ",
                  paste(names(df.split), collapse = ", ")))
    }
  } else if (size < 1) {
    n <- round(df.table * size, digits = 0)
  } else if (size >= 1) {
    if (all(df.table >= size) || isTRUE(replace)) {
      n <- setNames(rep(size, length.out = length(df.split)),
                    names(df.split))
    } else {
      message(
        "Some groups\n---",
        paste(names(df.table[df.table < size]), collapse = ", "),
        "---\ncontain fewer observations",
        " than desired number of samples.\n",
        "All observations have been returned from those groups.")
      n <- c(sapply(df.table[df.table >= size], function(x) x = size),
             df.table[df.table < size])
    }
  }
  temp <- lapply(
    names(df.split),
    function(x) df.split[[x]][sample(df.table[x],
                                     n[x], replace = F), ])
  set1 <- do.call("rbind", temp)
  set1
}

# k.fold description:
# 
# Takes a formula, dataset, number of folds and the output vector. 
# Uses nnet::multinom (see nnet v7.3-14 documentation - multinom).
# Divides data into folds, each time using a single fold as the test set
# and the rest as training. 
# Creates a predicted probalities and classes dataframe for each run 
# which is then combined to form a complete table.
# Returns a classifictaion table for predicted outcome vs. experimental observations,
# an accuracy measure based on that table and the probabilities table. 

#' Perform a k-fold cv for an ordinal logistic model  
#'
#' @param formula model formula
#' @param data a data frame with a class column 
#' @param folds number of folds (k)
#' @param outcome.column column with classes
#' @param stratify logical, should stratify in breaking to folds
#' @param sample.vector what is the ratio in sampling
#' between classes (for unbalanced sets) e.g. c(2,2,1)
#' @keywords internal
#' @return model stats
k.fold.log.ordinal <- function(formula, data, folds = NULL,
                               outcome.column = which(colnames(data) == 'class'),
                               stratify = F,
                               sample.vector = floor(round(summary(data$class)/min(summary(data$class))))) {
  models <- list()
  probalities <- list()
  acc <- list()
  class.pred <- list()
  if (stratify == T && !is.null(sample.vector) && is.null(folds)) {
    sets <- list()
    for (i in 1:ceiling(dim(data)[1] / sum(sample.vector))) {
      if (length(sets) == 0) {
        pool.data <- data
      } else {
        pool.data <- dplyr::setdiff(data, do.call(rbind, sets))
      }
      if (dim(pool.data)[1] < sum(sample.vector) & dim(pool.data)[1] > 0) {
        sets[[i]] <- pool.data
        pool.data <- data.frame()
      }
      if (nrow(pool.data) > 0 & any(summary(pool.data$class) < sample.vector)) {
        for (j in 1:nrow(pool.data)) {
          ifelse(j <= length(sets), sets[[j]] <- rbind(sets[[j]], pool.data[j, ]),
                 sets[[(j - length(sets))]] <- rbind(sets[[(j - length(sets))]], pool.data[j, ]))
        }
        is_empty <- function(x) (nrow(x) == 0 || is.null(x))
        sets <- sets[sapply(sets, is_empty) == FALSE]
      } else {
        if (nrow(pool.data) > 0) suppressWarnings(sets[[i]] <- stratified(pool.data, 'class', sample.vector))
      }
      if (i == ceiling(dim(data)[1] / sum(sample.vector))) {
        check <- dim(do.call(rbind, sets)) == dim(data)
        if (!all(check)) {
          print("Wasn't able to stritify as wanted. Check with another sample vector.")
        }
      }
    }
    for (i in 1:length(sets)) {
      sets[[i]] <- dplyr::mutate(sets[[i]], split.assign = i)
    }
    new_dat <- do.call(rbind, sets)
    new_dat <- dplyr::arrange(new_dat, flag)
    for (i in 1:length(sort(unique(new_dat$split.assign)))) {
      train <- new_dat[new_dat$split.assign != i,]
      test <- new_dat[new_dat$split.assign == i,]
      num.of.vars <- stringi::stri_count(formula, fixed = '+')
      start <- c(rep(0, (num.of.vars + 2)), 1)
      models[[match(i,sort(unique(new_dat$split.assign)))]] <- MASS::polr(formula,
                                                                          data = train,
                                                                          Hess = T, start = start,
                                                                          control = list(maxit = 1000))
      probalities[[match(i,sort(unique(new_dat$split.assign)))]] <- data.frame(predict(models[[match(i,sort(unique(new_dat$split.assign)))]], 
                                                                                       newdata = test,
                                                                                       "probs")*100)
      class.pred[[match(i,sort(unique(new_dat$split.assign)))]] <- data.frame(predict(models[[match(i,sort(unique(new_dat$split.assign)))]],
                                                                                      newdata = test,
                                                                                      "class"))
      probalities[[match(i,sort(unique(new_dat$split.assign)))]] <- cbind(probalities[[match(i,sort(unique(new_dat$split.assign)))]],
                                                                          class.pred[[match(i,sort(unique(new_dat$split.assign)))]],
                                                                          test$flag)
      names(probalities[[match(i,sort(unique(new_dat$split.assign)))]])[dim(probalities[[1]])[2] - 1] <- 'prediction'
      names(probalities[[match(i,sort(unique(new_dat$split.assign)))]])[dim(probalities[[1]])[2]] <- 'flag'
      
    }
  } else {
    if (folds == nrow(data)) {
      split.assign <- sample(1:folds, nrow(data), replace = F)
    } else {
      split.assign <- caret::createFolds(1:dim(data)[1], folds, list = F)
    }
    new_dat <- cbind(data, split.assign)
    for (i in 1:folds) {
      train <- new_dat[split.assign != i,]
      test <- new_dat[split.assign == i,]
      num.of.vars <- stringi::stri_count(formula, fixed = '+')
      start <- c(rep(0, (num.of.vars + 2)), 1)
      models[[match(i,1:folds)]] <- MASS::polr(formula,
                                               data = train,
                                               Hess = T, start = start,
                                               control = list(maxit = 1000))
      if (folds == nrow(data)) {
        probalities[[match(i,1:folds)]] <- data.table::transpose(data.frame(predict(models[[match(i,1:folds)]], 
                                                                                    newdata = test, "probs")*100))
      } else {
        probalities[[match(i,1:folds)]] <- data.frame(predict(models[[match(i,1:folds)]], 
                                                              newdata = test, "probs")*100)
      }
      class.pred[[match(i,1:folds)]] <- data.frame(predict(models[[match(i,1:folds)]], newdata = test, "class"))
      probalities[[match(i,1:folds)]] <- cbind(probalities[[match(i,1:folds)]],class.pred[[match(i,1:folds)]],test$flag)
      names(probalities[[match(i,1:folds)]])[dim(probalities[[1]])[2] - 1] <- 'prediction' 
      names(probalities[[match(i,1:folds)]])[dim(probalities[[1]])[2]] <- 'flag'
    }
  }
  
  probs <- data.frame(do.call(rbind, probalities))
  probs <- probs[order(probs$flag),]
  probs[,1:(dim(probalities[[1]])[2] - 2)] <- round(probs[,1:(dim(probalities[[1]])[2] - 2)],digits = 0)
  pred <- probs[,(dim(probalities[[1]])[2] - 1)]
  actual <- data[[outcome.column]]
  ct <- table(actual, pred)
  acc <- round((sum(diag(ct))/sum(ct))*100,2)
  ct.df <- data.frame(ct)
  TP <- ct.df$Freq[ct.df$actual == ct.df$pred]
  TP_FN <- list()
  TN <- list()
  TN_FP <- list()
  J <- list()
  for (i in levels(ct.df$actual)) {
    TP_FN[[i]] <- ct.df$Freq[ct.df$actual == i]
    TP_FN[[i]] <- sum(TP_FN[[i]])
    TP_FN[[i]] <- TP[which(levels(ct.df$actual) == i)]/TP_FN[[i]]
    TN[[i]] <- sum(ct.df$Freq[dplyr::intersect(which(ct.df$pred != i),
                                               which(ct.df$actual != i))])
    TN_FP[[i]] <- sum(ct.df$Freq[(ct.df$actual != i)])
    TN_FP[[i]] <- TN[[i]]/TN_FP[[i]]
    J[[i]] <- TP_FN[[i]] + TN_FP[[i]] - 1
  }
  J.small <- J[which.min(summary(data$class))]
  return(list(acc, J.small, ct, probs))
}

#' Perform a k-fold cv for a logistic model  
#'
#' @param formula model formula
#' @param data a data frame with a class column 
#' @param folds number of folds (k)
#' @param outcome.column column with classes
#' @param stratify logical, should stratify in breaking to folds
#' @param sample.vector what is the ratio in sampling
#' between classes (for unbalanced sets) e.g. c(2,2,1)
#' @keywords internal
#' @return model stats
k.fold.logistic <- function(formula, data, folds = NULL,
                            outcome.column = which(colnames(data) == 'class'),
                            stratify = F,
                            sample.vector = floor(round(summary(data$class)/min(summary(data$class))))) {
  models <- list()
  probalities <- list()
  acc <- list()
  class.pred <- list()
  if (stratify == T && !is.null(sample.vector) && is.null(folds)) {
    sets <- list()
    for (i in 1:ceiling(dim(data)[1] / sum(sample.vector))) {
      if (length(sets) == 0) {
        pool.data <- data
      } else {
        pool.data <- dplyr::setdiff(data, do.call(rbind, sets))
      }
      if (dim(pool.data)[1] < sum(sample.vector) & dim(pool.data)[1] > 0) {
        sets[[i]] <- pool.data
        pool.data <- data.frame()
      }
      if (nrow(pool.data) > 0 & any(summary(pool.data$class) < sample.vector)) {
        for (j in 1:nrow(pool.data)) {
          ifelse(j <= length(sets), sets[[j]] <- rbind(sets[[j]], pool.data[j, ]),
                 sets[[(j - length(sets))]] <- rbind(sets[[(j - length(sets))]], pool.data[j, ]))
        }
        is_empty <- function(x) (nrow(x) == 0 || is.null(x))
        sets <- sets[sapply(sets, is_empty) == FALSE]
      } else {
        if (nrow(pool.data) > 0) suppressWarnings(sets[[i]] <- stratified(pool.data, 'class', sample.vector))
      }
      if (i == ceiling(dim(data)[1] / sum(sample.vector))) {
        check <- dim(do.call(rbind, sets)) == dim(data)
        if (!all(check)) {
          print("Wasn't able to stritify as wanted. Check with another sample vector.")
        }
      }
    }
    for (i in 1:length(sets)) {
      sets[[i]] <- dplyr::mutate(sets[[i]], split.assign = i)
    }
    new_dat <- do.call(rbind, sets)
    new_dat <- dplyr::arrange(new_dat, flag)
    for (i in 1:length(sort(unique(new_dat$split.assign)))) {
      train <- new_dat[new_dat$split.assign != i,]
      test <- new_dat[new_dat$split.assign == i,]
      num.of.vars <- stringi::stri_count(formula, fixed = '+')
      start <- c(rep(0, (num.of.vars + 2)), 1)
      models[[match(i,sort(unique(new_dat$split.assign)))]] <- nnet::multinom(formula,
                                                                              data = train,
                                                                              maxit = 2000,
                                                                              trace = FALSE)
      probalities[[match(i,sort(unique(new_dat$split.assign)))]] <- data.frame(predict(models[[match(i,sort(unique(new_dat$split.assign)))]], 
                                                                                       newdata = test,
                                                                                       "probs")*100)
      class.pred[[match(i,sort(unique(new_dat$split.assign)))]] <- data.frame(predict(models[[match(i,sort(unique(new_dat$split.assign)))]],
                                                                                      newdata = test,
                                                                                      "class"))
      probalities[[match(i,sort(unique(new_dat$split.assign)))]] <- cbind(probalities[[match(i,sort(unique(new_dat$split.assign)))]],
                                                                          class.pred[[match(i,sort(unique(new_dat$split.assign)))]],
                                                                          test$flag)
      names(probalities[[match(i,sort(unique(new_dat$split.assign)))]])[dim(probalities[[1]])[2] - 1] <- 'prediction'
      names(probalities[[match(i,sort(unique(new_dat$split.assign)))]])[dim(probalities[[1]])[2]] <- 'flag'
      
    }
  } else {
    if (folds == nrow(data)) {
      split.assign <- sample(1:folds, nrow(data), replace = F)
    } else {
      split.assign <- caret::createFolds(1:dim(data)[1], folds, list = F)
    }
    new_dat <- cbind(data, split.assign)
    for (i in 1:folds) {
      train <- new_dat[split.assign != i,]
      test <- new_dat[split.assign == i,]
      num.of.vars <- stringi::stri_count(formula, fixed = '+')
      start <- c(rep(0, (num.of.vars + 2)), 1)
      models[[match(i,1:folds)]] <- nnet::multinom(formula,
                                                   data = train,
                                                   maxit = 2000,
                                                   trace = FALSE)
      if (folds == nrow(data)) {
        probalities[[match(i,1:folds)]] <- data.table::transpose(data.frame(predict(models[[match(i,1:folds)]], 
                                                                                    newdata = test, "probs")*100))
      } else {
        probalities[[match(i,1:folds)]] <- data.frame(predict(models[[match(i,1:folds)]], 
                                                              newdata = test, "probs")*100)
      }
      class.pred[[match(i,1:folds)]] <- data.frame(predict(models[[match(i,1:folds)]], newdata = test, "class"))
      probalities[[match(i,1:folds)]] <- cbind(probalities[[match(i,1:folds)]],class.pred[[match(i,1:folds)]],test$flag)
      names(probalities[[match(i,1:folds)]])[dim(probalities[[1]])[2] - 1] <- 'prediction' 
      names(probalities[[match(i,1:folds)]])[dim(probalities[[1]])[2]] <- 'flag'
    }
  }
  
  probs <- data.frame(do.call(rbind, probalities))
  probs <- probs[order(probs$flag),]
  probs[,1:(dim(probalities[[1]])[2] - 2)] <- round(probs[,1:(dim(probalities[[1]])[2] - 2)],digits = 0)
  pred <- probs[,(dim(probalities[[1]])[2] - 1)]
  actual <- data[[outcome.column]]
  ct <- table(actual, pred)
  acc <- round((sum(diag(ct))/sum(ct))*100,2)
  ct.df <- data.frame(ct)
  TP <- ct.df$Freq[ct.df$actual == ct.df$pred]
  TP_FN <- list()
  TN <- list()
  TN_FP <- list()
  J <- list()
  for (i in levels(ct.df$actual)) {
    TP_FN[[i]] <- ct.df$Freq[ct.df$actual == i]
    TP_FN[[i]] <- sum(TP_FN[[i]])
    TP_FN[[i]] <- TP[which(levels(ct.df$actual) == i)]/TP_FN[[i]]
    TN[[i]] <- sum(ct.df$Freq[dplyr::intersect(which(ct.df$pred != i),
                                               which(ct.df$actual != i))])
    TN_FP[[i]] <- sum(ct.df$Freq[(ct.df$actual != i)])
    TN_FP[[i]] <- TN[[i]]/TN_FP[[i]]
    J[[i]] <- TP_FN[[i]] + TN_FP[[i]] - 1
  }
  J.small <- J[which.min(summary(data$class))]
  return(list(acc, J.small, ct, probs))
}


#' Get ordinal logistic models of all feature subsets within a scope  
#'
#' @param data a data frame with a class column 
#' @param out.col what column holds the class
#' @param min minimum number of variables
#' @param max maximum number of variables
#' @return list of models, ranked by McFadden's pseudo R squared
#' @export
model.subset.log.ordinal <- function(data, out.col = which(colnames(data) == 'class'),
                              min = 1, max = floor(dim(data)[1]/5)) {
  output <- stringr::str_c('`',names(data[out.col]),'`')
  vars <- names(data[,-out.col])
  vars <- vars[vars != 'flag']
  for (i in 1:length(vars)) {
    vars[i] <- stringr::str_c('`',vars[i],'`')
  }
  comb.list <- list()
  R2.list <- list()
  for (i in min:max) {
    comb.list[[i]] <- data.frame(aperm(combn(vars,i)), stringsAsFactors = F)
    comb.list[[i]][, dim(comb.list[[i]])[2] + 1] <- do.call(paste, 
                                                            c(comb.list[[i]][names(comb.list[[i]])],
                                                              sep = " + "))
    names(comb.list[[i]])[dim(comb.list[[i]])[2]] <- 'formula'
    for (co in names(comb.list[[i]])[1:length(names(comb.list[[i]])) - 1]) comb.list[[i]][co] <- NULL
  }
  comb.list <- plyr::compact(comb.list)
  forms <- do.call(rbind,comb.list)
  names(forms) <- 'formula'
  for (i in 1:dim(forms)[1]) {
    forms$formula[i] <- stringr::str_c(output,' ~ ',forms$formula[i])
    num.of.vars <- stringi::stri_count(forms$formula[i], fixed = '+')
    start <- c(rep(0, num.of.vars + 2), 1)
    test.1 <- MASS::polr(forms$formula[i],
                         data = data,
                         Hess = T, start = start,
                         control = list(maxit = 100))
    test.0 <- MASS::polr(class ~ 1,data = data, Hess = T)
    loglik.0 <- -test.0$deviance/2
    loglik.1 <- -test.1$deviance/2
    R2.list[i] <- round((1 - (loglik.1/loglik.0)),digits = 3)
  }
  forms[,2] <- do.call(rbind,R2.list)
  names(forms)[2] <- 'McFadden R2'
  out.models <- head(dplyr::arrange(forms,desc(forms[,2])),10)
  return(out.models)
}

#' Get logistic models of all feature subsets within a scope  
#'
#' @param data a data frame with a class column 
#' @param out.col what column holds the class
#' @param min minimum number of variables
#' @param max maximum number of variables
#' @return list of models, ranked by McFadden's pseudo R squared
#' @export
model.subset.logistic <- function(data, out.col = which(colnames(data) == 'class'),
                                  min = 1, max = floor(dim(data)[1]/5)) {
  output <- stringr::str_c('`',names(data[out.col]),'`')
  vars <- names(data[,-out.col])
  vars <- vars[vars != 'flag']
  for (i in 1:length(vars)) {
    vars[i] <- stringr::str_c('`',vars[i],'`')
  }
  comb.list <- list()
  R2.list <- list()
  for (i in min:max) {
    comb.list[[i]] <- data.frame(aperm(combn(vars,i)), stringsAsFactors = F)
    comb.list[[i]][, dim(comb.list[[i]])[2] + 1] <- do.call(paste, 
                                                            c(comb.list[[i]][names(comb.list[[i]])],
                                                              sep = " + "))
    names(comb.list[[i]])[dim(comb.list[[i]])[2]] <- 'formula'
    for (co in names(comb.list[[i]])[1:length(names(comb.list[[i]])) - 1]) comb.list[[i]][co] <- NULL
  }
  comb.list <- plyr::compact(comb.list)
  forms <- do.call(rbind,comb.list)
  names(forms) <- 'formula'
  for (i in 1:dim(forms)[1]) {
    forms$formula[i] <- stringr::str_c(output,' ~ ',forms$formula[i])
    test.1 <- nnet::multinom(forms$formula[i],
                             data = data,
                             maxit = 2000,
                             trace = FALSE)
    test.0 <- nnet::multinom(class ~ 1,
                             data = data,
                             maxit = 2000,
                             trace = FALSE)
    loglik.0 <- -test.0$deviance/2
    loglik.1 <- -test.1$deviance/2
    R2.list[i] <- round((1 - (loglik.1/loglik.0)),digits = 3)
  }
  forms[,2] <- do.call(rbind,R2.list)
  names(forms)[2] <- 'McFadden R2'
  out.models <- head(dplyr::arrange(forms,desc(forms[,2])),10)
  return(out.models)
}

# kf.iter description:
# 
# Takes a formula, dataset, number of folds and the output vector along with 
# number of iterations. 
# Runs k.fold.ordinal.log by iterations numbers and returns the average 
# accuracy over all iterations, the best and worst accuracy and 
# best and worst classification tables.

#' Iterate a k-fold cv - ordinal model 
#'
#' @param formula model formula
#' @param data a data frame with a class column 
#' @param folds number of folds (k)
#' @param out.col column with classes
#' @param stratify logical, should stratify in breaking to folds
#' @param sample.vector what is the ratio in sampling
#' between classes (for unbalanced sets) e.g. c(2,2,1)
#' @param iterations number of iterations
#' @param verbose should spit out information 
#' @return model stats
#' @keywords internal
kf.iter.log.ordinal <- function(formula, data, folds = NULL, out.col=which(colnames(data) == 'class'),
                            stratify = F,
                            sample.vector = floor(round(summary(data$class)/min(summary(data$class)))),
                            iterations, verbose = F) {
  iter.list <- list()
  ct.list <- list()
  for (i in 1:iterations) {
    mod <- k.fold.log.ordinal(formula, data, folds, out.col, stratify, sample.vector)
    iter.list[[match(i,1:iterations)]] <- mod[[1]]
    ct.list[[match(i,1:iterations)]] <- mod[[3]]
  }
  over.all.accuracy <- round(Reduce(`+`,iter.list)/iterations,digits = 2)
  best <- iter.list[which.max(iter.list)]
  worst <- iter.list[which.min(iter.list)]
  Accuracies <- knitr::kable(cbind(over.all.accuracy, best,worst))
  tab <- suppressMessages(cbind(reshape2::dcast(data.frame(ct.list[which.max(iter.list)]),actual~pred),
                                rep("***",3),
                                reshape2::dcast(data.frame(ct.list[which.min(iter.list)]),actual~pred)))
  names(tab)[5] <- ''
  cts <- knitr::kable(tab)
  if (verbose == T) { 
    print(cts)
    print(Accuracies)
  }
  invisible(over.all.accuracy)
}

#' Iterate a k-fold cv - logistic model 
#'
#' @param formula model formula
#' @param data a data frame with a class column 
#' @param folds number of folds (k)
#' @param out.col column with classes
#' @param stratify logical, should stratify in breaking to folds
#' @param sample.vector what is the ratio in sampling
#' between classes (for unbalanced sets) e.g. c(2,2,1)
#' @param iterations number of iterations
#' @param verbose should spit out information 
#' @return model stats
#' @keywords internal
kf.iter.logistic <- function(formula, data, folds = NULL, out.col=which(colnames(data) == 'class'),
                             stratify = F,
                             sample.vector = floor(round(summary(data$class)/min(summary(data$class)))),
                             iterations, verbose = F) {
  iter.list <- list()
  ct.list <- list()
  for (i in 1:iterations) {
    mod <- k.fold.logistic(formula, data, folds, out.col, stratify, sample.vector)
    iter.list[[match(i,1:iterations)]] <- mod[[1]]
    ct.list[[match(i,1:iterations)]] <- mod[[3]]
  }
  over.all.accuracy <- round(Reduce(`+`,iter.list)/iterations,digits = 2)
  best <- iter.list[which.max(iter.list)]
  worst <- iter.list[which.min(iter.list)]
  Accuracies <- knitr::kable(cbind(over.all.accuracy, best,worst))
  tab <- suppressMessages(cbind(reshape2::dcast(data.frame(ct.list[which.max(iter.list)]),actual~pred),
                                rep("***",3),
                                reshape2::dcast(data.frame(ct.list[which.min(iter.list)]),actual~pred)))
  names(tab)[5] <- ''
  cts <- knitr::kable(tab)
  if (verbose == T) { 
    print(cts)
    print(Accuracies)
  }
  invisible(over.all.accuracy)
}


# mod.info description:
# 
# Takes a logistic model object.
# using the full data set as both training and testing sets, model information is extracted.
# Creates a class.table object to be used in confusion matrix plotting,
# computes and returns variable importance as the sum of each coefficient's absolute value
# over all predicted classes (See caret v6.0-86 documentation - varImp),
# accuracy based on classification table, 
# coefficient's exponent as another measure of importance 
# and McFadden's pseudo R squared for multinomial logistic regression (T. Domencich and D. L. McFadden, 1975).

#' Model summary - ordinal model 
#'
#' @param model model object - MASS::polr object
#' @param data data frame with class column 
#' @return model stats
#' @export
mod.info.log.ordinal <- function(model,data) {
  pred <- predict(model, newdata = data, 'class')
  actual <- data$class
  class.table <- table(actual, pred)
  
  Accuracy <- paste(round((sum(diag(class.table))/sum(class.table))*100,2),"%",sep = '')
  test.0 <- MASS::polr(class ~ 1,data = data, Hess = T)
  test.1 <- model
  
  loglik.0 <- -test.0$deviance/2
  loglik.1 <- -test.1$deviance/2
  
  McFadden_R2 <- round((1 - (loglik.1/loglik.0)),digits = 3)
  st <- cbind(Accuracy,McFadden_R2)
  
  stats <- knitr::kable(st)
  print(stats)
  
  sum.mod <- summary(model)
  coefs.mod <- knitr::kable(sum.mod$coefficients)
  print(coefs.mod)
  return(class.table)
}

#' Model summary - logistic model 
#'
#' @param model model object - nnet::multinom object
#' @param data data frame with class column 
#' @return model stats
#' @export
mod.info.logistic <- function(model, data) {
  pred <- predict(model, newdata = data, 'class')
  actual <- data$class
  class.table <- table(actual, pred)
  
  Accuracy <- paste(round((sum(diag(class.table))/sum(class.table))*100,2),"%",sep = '')
  test.0 <- nnet::multinom(class ~ 1,
                           data = data,
                           maxit = 2000,
                           trace = FALSE)
  test.1 <- model
  
  loglik.0 <- -test.0$deviance/2
  loglik.1 <- -test.1$deviance/2
  
  McFadden_R2 <- round((1 - (loglik.1/loglik.0)),digits = 3)
  st <- cbind(Accuracy,McFadden_R2)
  
  stats <- knitr::kable(st)
  print(stats)
  
  test.1$call$formula <- as.formula(paste('class',
                                          '~',
                                          paste(test.1$coefnames[-1],
                                                collapse = '+')))
  test.1$call$data <- data
  sum.mod <- summary(test.1)
  coefs.mod <- knitr::kable(sum.mod$coefficients)
  print(coefs.mod)
  return(class.table)
}

# ct.plot description:
# 
# Takes a class table and computes;
# * Recall - TP/(TP+FN) & complimentary percentages of wrong prdicitions in each cell.
# * Precision - TP/(TP+FP) for each prediction class.
# * Accuracy - all.TP/total.
# Classifies into Right and Wrong predictions and presents class size.

#' Create a nice classification table (confusion matrix)
#' 
#' Takes a mod.info object (x <- mod.info(blah) and then 
#' ct.plot(x))
#'
#' @param class.table a mod.info result
#' @return classification table figure
#' @export
ct.plot <- function(class.table) {
  ct <- as.matrix(class.table)
  classes <- nrow(ct)
  total <- rep(NA, nrow(ct))
  ct <- cbind(ct,total)
  for (i in 1:dim(ct)[1]) {
    ct[i, (dim(ct)[1] + 1)] <- sum(ct[i, 1:(dim(ct)[2] - 1)])
  }
  total <- rep(NA, ncol(ct))
  ct <- rbind(ct,total)
  for (i in 1:dim(ct)[2]) {
    ct[dim(ct)[2], i] <- sum(ct[1:(dim(ct)[1] - 1), i])
  }
  ct <- reshape2::melt(ct)
  names(ct) <- c('Exp','Pred','Freq')
  ct$Exp <- as.factor(ct$Exp)
  ct$Pred <- as.factor(ct$Pred)
  
  ct <- dplyr::mutate(ct, Right.Wrong=rep(NA, nrow(ct)))
  ct <- dplyr::mutate(ct,prop = rep(NA, nrow(ct)))
  
  
  
  for (i in as.numeric(row.names(ct[ct$Exp != 'total' & ct$Pred != 'total', ]))) {
    for (j in 1:classes) {
      if (ct$Exp[[i]] == as.character(j)) {
        ct$prop[[i]] <- ct$Freq[[i]]/sum(ct$Freq[ct$Exp==as.character(j) & ct$Pred!='total'])
      }
    }
    if (ct$Exp[i] == ct$Pred[i]) {
      ct$Right.Wrong[i] <- 'True'
    } else {
      ct$Right.Wrong[i] <- 'False'
    }
  } 
  ct$prop <- round(ct$prop,digits = 3)*100
  
  for (i in as.numeric(row.names(ct[ct$Exp =='total' | ct$Pred =='total', ]))) {
    if (ct$Pred[i]=='total' & ct$Exp[i]!='total') {
      ct$prop[i] <- (ct$Freq[i]/ct$Freq[ct$Exp == 'total' & ct$Pred == 'total'])*100
      ct$Right.Wrong[i] <- 'Size'
    }
    if (ct$Pred[i]!='total' & ct$Exp[i]=='total') {
      ct$prop[[i]] <- (ct$Freq[ct$Exp == ct$Pred[i] & ct$Pred == ct$Pred[i]]/ct$Freq[i])*100
      ct$Right.Wrong[i] <- 'Precision'
    }
    if (ct$Pred[i]=='total' & ct$Exp[i]=='total') {
      ct$prop[[i]] <- (sum(diag(class.table))/ct$Freq[i])*100
      ct$Right.Wrong[i] <- 'Accuracy'
    }
  } 
  ct$prop <- round(ct$prop,digits = 1)
  
  ct <- dplyr::mutate(ct,value=ct$prop)
  
  for (i in 1:nrow(ct)) {
    if (ct$Right.Wrong[i]!='True' & ct$Right.Wrong[i]!='False') {
      ct$prop[i] <- 45
    }
  }
  
  base <- ggplot2::ggplot(data = ct, 
                          mapping = ggplot2::aes(x =ordered(Pred, levels =sort(unique(Pred))),
                                                 y=ordered(Exp, levels =rev(sort(unique(Exp)))),
                                                 fill = Right.Wrong,
                                                 alpha = prop))+
    ggplot2::geom_tile(color = 'black',size = 1.2)+
    ggplot2::coord_fixed()+
    ggplot2::geom_text(ggplot2::aes(label=paste(Freq,"\n",
                                                '(',value,'%',')',sep = '')),
                       size=5,vjust = .5, fontface  = "bold", alpha = 1) +
    ggplot2::scale_fill_manual(values = c(True = '#62a188', False = '#d1505b',
                                          Size = 'lightgrey', Precision = '#FFCC33',
                                          Accuracy = 'lightblue'))+
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 11,face = 'bold'),
                   axis.text.y = ggplot2::element_text(size = 11,face = 'bold'),
                   axis.title.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())
  base
}

# prob.heatmap description:
# 
# Takes a nnet::multinom model object.
# Represents the final probabilities table as a heat map with a color gradient based
# on probability percentage. The correctly predicted case names are colored green,
# wrong but second in percentage are colored orange and absolute wrong is colored red.
# The experimental outcome is presented in its own column, color coded for convenience.

#' Create a nice probabilities of predictions table
#'
#' @param model a classifier model (nnet/MASS)
#' @param data a data frame with a class column
#' @return probabilities table figure
#' @export
prob.heatmap <- function(model,data) {
  classes <- length(unique(data$class))
  pred <- predict(model,newdata = data, 'class')
  r.w <- pred == data$class
  probs <- fitted(model)*100
  verif <- data.frame(cbind(data$class, pred, r.w, probs, rep(NA, nrow(probs))))
  row.names(verif) <- row.names(probs)
  colnames(verif)[c(1, dim(verif)[2])] <- c('Exp','color')
  for (i in 1:dim(verif)[1]) {
    second <- sort(as.numeric(verif[i, 4:(4 + classes)]), decreasing = T)[2]
    where.is.sec <- which(as.numeric(verif[i, 4:(4 + (classes - 1))]) == second)
    if (verif$r.w[i] == 1) { 
      verif$color[i] <- "#66a180" # green
    } else {
      if (as.numeric(verif$Exp[i]) == where.is.sec) {
        verif$color[i] <- 'tan1'
      } else {
        verif$color[i] <- '#d1505b'
      }
    }
  }
  
  pro.df <- data.frame(probs)
  pro.df <- tibble::rownames_to_column(pro.df)
  pro.df[, 1] <- factor(pro.df[, 1], levels = pro.df[, 1])
  pro.df[, (dim(pro.df)[2] + 1)] <- as.numeric(data$class)
  colnames(pro.df) <- c('Name',as.character(1:classes),'Exp')
  row.names(pro.df) <- row.names(probs)
  
  long <- reshape2::melt(pro.df, id.vars = 'Name')
  long[,3] <- round(long[,3],digits = 2)
  long <- dplyr::mutate(long, exp_shape = rep(NA,nrow(long)))
  
  
  for (i in 1:nrow(long)) {
    for (j in 1:classes) {
      if (long$variable[i] == 'Exp' & long$value[i] == j) {
        long$exp_shape[i] <- j
      }
    }
  }
  col_vec <- vector(mode = 'numeric')
  coloring <- c("darkgoldenrod4", "slateblue", 'darksalmon',
                                'darkblue', 'navajowhite4',
                                'darkcyan', 'chocolate4',
                                "coral3","cornsilk4",'darkslateblue')
                                for (i in 1:length(long[long$variable == 'Exp',4])) {
                                  col_vec[i] <- coloring[long[long$variable == 'Exp',4][i]]
                                }
  
  label.vec <- as.character(long[long$variable == 'Exp',4])
  
  prob.heatmap <- ggplot2::ggplot(mapping = ggplot2::aes(x = variable,
                                                         y = ordered(Name, 
                                                                    levels = rev(factor(pro.df$Name, 
                                                                                       levels = pro.df$Name))))) +
    ggplot2::geom_tile(data = long[long$variable != 'Exp',], 
                       color = 'black', ggplot2::aes(fill = value))+
    ggplot2::coord_fixed(ratio = 0.2)+
    ggplot2::geom_text(data = long[long$variable != 'Exp',], 
                       ggplot2::aes(label=value))+
    ggplot2::scale_fill_gradient(name = "% probability",
                                 low = "#FFFFFF",
                                 high = "dodgerblue2",
                                 guide = ggplot2::guide_colorbar(frame.colour = "black", 
                                                                 ticks.colour = "black"))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, face = 'bold',),
                   axis.text.y = ggplot2::element_text(size = 10, face = 'bold', 
                                                       colour = rev(verif$color)),
                   axis.title.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())+
    ggplot2::scale_x_discrete(position = "top",limits = levels(long$variable))+
    ggplot2::geom_tile(data = long[long$variable == 'Exp',],
                       alpha = 0, inherit.aes = F,
                       ggplot2::aes(x=rev(variable),
                                    y=ordered(Name, levels =rev(factor(pro.df$Name, 
                                                                       levels = pro.df$Name)))))+
    ggplot2::geom_text(data = long[long$variable == 'Exp',], label = label.vec,
                       size = 4, color = col_vec, fontface = 'bold')
  
  prob.heatmap
}

# simi.sampler description:
#
# Takes data and class as input, yields a subset of the class based
# on a cosine similarity test of the class samples with the class mean vector. 
# This tests how "similar" each sample is to the class itself. 
# The sampler takes the smallest class as the size of subset and measures the 
# the distances of all samples in the class with "poles" stationed in equal spacing 
# from the most similar to the least. 

#' a similarity based sampling algorithm 
#' 
#' Used to resize (under-sample) classification oriented data sets.
#' Returns the samples row numbers.  
#'
#' @param data data frame with class column
#' @param class class of interest number
#' @param compare.with choose the class with which similarity is computed.
#' Defaults to 0, being similarity of each sample with its own group. 
#' Any other number will compare with the class represented by that number. 
#' @param plot create a plot of the similarity before and after sampling.
#' @param sample.size how many are to be sampled. 
#' Defaults to the number of sample in the smallest class. 
#' @return A vector with row numbers of sampled observations
#' @export
simi.sampler <- function(data, class,
                         compare.with = 0,
                         plot = F,
                         sample.size = min(summary(as.factor(data$class)))) {
  out.col <- which(colnames(data) == 'class')
  vars <- names(data[,-out.col])
  vars <- vars[vars != 'flag']
  vars <- vars[vars != 'tag']
  sampler.data <- data[vars]
  sampler.data <- data.frame(apply(sampler.data, 2, as.numeric))
  sampler.data <- data.frame(scale(sampler.data, T, T))
  var.nums <- 1:dim(sampler.data)[2]
  classes <- length(unique(data$class))
  
  # compute mean values vector for each of the groups - for similarity
  for (i in 1:classes) {
    assign(paste('class.', as.character(i), '.vector', sep = ''), 
           apply(sampler.data[data$class == i,][,vars], 2, mean))
    assign(paste('class.', as.character(i), '.mag', sep = ''), 
           sqrt(sum(apply(sampler.data[data$class == i,][,vars], 2, mean)^2)))
  }
  
  # compute similarity of each group with itself - same class
  new.col <- dim(sampler.data)[2] + 1
  for (r in 1:nrow(sampler.data)) {
    for (i in 1:classes) {
      if (data$class[r] == i) {
        vec <- get(ls()[which(grepl(paste('.',as.character(i),'.vector',sep = ''), ls()))])
        mag <- get(ls()[which(grepl(paste('.',as.character(i),'.mag',sep = ''), ls()))]) 
        sampler.data[r ,new.col] <- sum(vec*sampler.data[r, 1:(new.col - 1)])/((mag*sqrt(sum(sampler.data[r, 1:(new.col - 1)]^2))))
      }
    }
  }
  
  # compute similarity between classes
  simi.df <- data.frame(matrix(ncol = classes, nrow = nrow(data)))
  for (i in 1:nrow(simi.df)) {
    for (j in 1:classes) {
      vec <- get(ls()[which(grepl(paste('.',as.character(j),'.vector',sep = ''), ls()))])
      mag <- get(ls()[which(grepl(paste('.',as.character(j),'.mag',sep = ''), ls()))]) 
      simi.df[i ,j] <- sum(vec*sampler.data[i, 1:(new.col - 1)])/((mag*sqrt(sum(sampler.data[i, 1:(new.col - 1)]^2))))
    }
  }
  names(simi.df) <- as.character(seq(1, classes))
  
  names(sampler.data)[new.col] <- 'class.similarity'
  simi.table <- cbind(sampler.data[new.col], simi.df)
  simi.table <- dplyr::mutate(simi.table, data$class)
  simi.table <- dplyr::mutate(simi.table, row.names(data))
  simi.table <- dplyr::mutate(simi.table, data$flag)
  colnames(simi.table)[c(1, (2 + classes):dim(simi.table)[2])] <- c('same.class','class','Name', 'flag')
  
  
  plot.sim.before <- ggplot2::ggplot(simi.table, 
                                     ggplot2::aes(simi.table[,compare.with + 1],
                                                  y = class, label = Name)) +
    ggplot2::geom_point(ggplot2::aes(color = class)) +
    ggrepel::geom_text_repel(ggplot2::aes(label = Name), max.overlaps = 25) +
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
    ggplot2::xlab(names(simi.table)[compare.with + 1]) + 
    ggplot2::ggtitle('Similarity Before Group Truncation')
  
  simi.class <- simi.table[simi.table$class == class, compare.with + 1]
  steps <-  sort(seq(min(simi.class), max(simi.class), (max(simi.class) - min(simi.class))/(sample.size - 1)))
  dis.mat <- matrix(ncol = length(steps), nrow = length(sort(simi.class)))
  for (i in 1:nrow(dis.mat)) {
    dis.mat[i, ] <- abs(simi.class[i] - steps)
  }
  pts <- vector()
  row.names(dis.mat) <- as.character(simi.table$flag[simi.table$class == class])
  if (length(steps) < length(simi.class)) {
    for (i in 1:ncol(dis.mat)) {
      drop <- which.min(dis.mat[, i])
      pts[i] <- simi.table[as.numeric(names(drop)), 1]
      dis.mat <- dis.mat[-drop, ]
    }
  } else {
    pts <- simi.table[as.numeric(row.names(dis.mat)), 1]
  }
  keep <- as.numeric(simi.table$flag[simi.table$same.class %in% pts])
  
  class.rows <- simi.table$flag[simi.table$class == class]
  if (min(class.rows) != 1) {
    truncated <- unique(c(1:(min(class.rows) + 1), keep, (max(class.rows) + 1):nrow(data)))
  } else {
    truncated <- unique(c(keep, (max(class.rows) + 1):nrow(data)))
  }
  simi.plot.data <- simi.table[truncated,][complete.cases(simi.table[truncated,]),]
  plot.sim.after <- ggplot2::ggplot(simi.plot.data,
                                    ggplot2::aes(x = simi.plot.data[, compare.with + 1],
                                                 y = class, label = Name)) +
    ggplot2::geom_point(ggplot2::aes(color = class)) +
    ggrepel::geom_text_repel(ggplot2::aes(label = Name), max.overlaps = 25) +
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
    ggplot2::xlab(names(simi.table)[compare.with + 1]) + 
    ggplot2::ggtitle('Similarity After Group Truncation')
  
  if (isTRUE(plot)) gridExtra::grid.arrange(plot.sim.before, plot.sim.after, ncol = 2)
  
  return(keep)
}


#' Generate a model report - ordinal logistic regression
#' 
#' A full, detailed report, includes plots  
#'
#' @param dataset data frame with class column
#' @param min minimum number of variables for model sub-setting 
#' @param max maximum number of variables for model sub-setting
#' @param model.number which model (out of the list) should be summarized
#' @return A full report
#' @export
model.report.log.ordinal <- function(dataset, min = 1,
                                     max = (min + 1), model.number = 1) {
  
  cat(tools::file_path_sans_ext(basename(dataset)))
  data <- data.frame(data.table::fread(dataset, header = T,
                                       check.names = T))
  RN <- data[,1]
  data <- data[,-1]
  data <- data[complete.cases(data), ]
  CN <- names(data)
  names(data) <- CN
  if (!'class' %in% CN) stop('One of the columns must be named class')
  row.names(data) <- RN
  data$class <- as.factor(data$class)
  
  models <-  model.subset.log.ordinal(data, min = min, max = max)
  
  knitr::kable(models)
  
  print(knitr::kable(models))
  
  test.form <- models[model.number, 1]
  
  test.data <- data
  
  num.of.vars <- stringi::stri_count(test.form, fixed = '+')
  start <- c(rep(0, num.of.vars + 2), 1)
  test <- MASS::polr(test.form,
                     data = test.data,
                     Hess = T, start = start)
  # LOOCV
  
  cat('
  Leave-one-out Cross Validation - single iteration
      ')
  kf.iter.log.ordinal(test.form, test.data, folds = nrow(test.data), iterations = 1, verbose = T)
  
  # size of smallest class-fold stratified CV
  cat('
  Smallest-group-fold Cross Validation (100 iterations)
  
  Classification tables for best (left) and worst (right) iterations. 
      ')
  kf.iter.log.ordinal(test.form, test.data, stratify = T, iterations = 100, verbose = T)
  
  cat('
  Model Summary
      ')
  ct <- mod.info.log.ordinal(test, test.data)
  
  CT <- ct.plot(ct)
  
  HM <- prob.heatmap(test, data)
  
  gridExtra::grid.arrange(CT, HM, ncol = 2)
}

#' Generate a model report - logistic regression
#' 
#' A full, detailed report, includes plots  
#'
#' @param dataset data frame with class column
#' @param min minimum number of variables for model sub-setting 
#' @param max maximum number of variables for model sub-setting
#' @param model.number which model (out of the list) should be summarized
#' @return A full report
#' @export
model.report.logistic <- function(dataset, min = 1,
                                  max = (min + 1), model.number = 1) {
  cat(tools::file_path_sans_ext(basename(dataset)))
  data <- data.frame(data.table::fread(dataset, header = T,
                                           check.names = T))
  RN <- data[,1]
  data <- data[,-1]
  data <- data[complete.cases(data), ]
  CN <- names(data)
  names(data) <- CN
  if (!'class' %in% CN) stop('One of the columns must be named class')
  row.names(data) <- RN
  data$class <- as.factor(data$class)
  
  models <-  model.subset.logistic(data, min = min, max = max)
  
  print(knitr::kable(models))
  
  test.form <- models[model.number, 1]
  
  test.data <- data
  
  test <- nnet::multinom(test.form,
                         data = test.data,
                         maxit = 2000, 
                         trace = F)
  # LOOCV
  cat('
  Leave-one-out Cross Validation - single iteration
      ')
  kf.iter.logistic(test.form, test.data, folds = nrow(test.data), iterations = 1, verbose = T)
  
  # size of smallest class-fold stratified CV
  cat('
  Smallest-group-fold Cross Validation (100 iterations)
  
  Classification tables for best (left) and worst (right) iterations. 
      ')
  kf.iter.logistic(test.form, test.data, stratify = T, iterations = 100, verbose = T)
  
  cat('
  Model Summary
      ')
  ct <- mod.info.logistic(test, test.data)
  
  CT <- ct.plot(ct)
  
  HM <- prob.heatmap(test, data)
  
  gridExtra::grid.arrange(CT, HM, ncol = 2)
}

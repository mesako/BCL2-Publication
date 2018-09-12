library(igraph)
library(flowCore)
library(robustbase)
library(ggplot2)
library(ggfortify)
library(GGally)
library(caret)
library(reshape2)
library(pROC)
library(glmnet)
library(randomForest)
library(glinternet)
library(preprocessCore)
library(FLOWMAPR)

CollapseBy <- function(named.list.of.df, labels.in.order,
                       new.col.name, subsample = FALSE) {
  single.df <- c()
  if (length(labels.in.order) != length(named.list.of.df)) {
    stop("Not enough labels for each df in list!")
  }
  for (i in 1:length(labels.in.order)) {
    this.df <- which(grepl(names(named.list.of.df), pattern = labels.in.order[i]))
    temp <- named.list.of.df[[this.df]]
    if (subsample != FALSE) {
      temp <- temp[sample(nrow(temp), subsample), ]
    }
    temp <- cbind(temp, rep(labels.in.order[i], times = nrow(temp)))
    if (is.factor(temp[, ncol(temp)])) {
      print("Converting from factor to character")
      temp[, ncol(temp)] <- as.character(temp[, ncol(temp)])
    }
    colnames(temp)[ncol(temp)] <- new.col.name
    single.df <- rbind(single.df, temp)
    rm(temp)
  }
  return(single.df)
}

ConvertClassto01 <- function(data, column, true.label) {
  if (length(grep(column, colnames(data))) == 0) {
    stop("column is not one of the columns in data!")
  }
  fixed.data <- data
  fixed.data[, column] <- fixed.data[, column] == true.label
  colnames(fixed.data)[grep(column, colnames(data))] <- true.label
  return(fixed.data)
}

SaveROC <- function(model, data, response, output.name,
                    ci, smooth) {
  prediction <- predict(model, newdata = data,
                        s = model$bestTune$lambda, type = "prob")
  prediction <- prediction[, 2]
  this.roc <- roc(data[, response], prediction, percent = TRUE, algorithm = 1,
                  quiet = TRUE, smooth = smooth, auc = TRUE, ci = ci, plot = TRUE)
  saveRDS(this.roc, file = output.name)
  base.name <- unlist(strsplit(output.name, split = "\\."))[1]
  png(paste(base.name, ".png", sep = ""), width = 4,
      height = 4, units = "in", res = 400)
  plot(this.roc, print.auc = TRUE, col = "red", print.auc.y = 5,
       print.auc.x = 85)
  dev.off()
  return(this.roc)
}

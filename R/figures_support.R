library(igraph)
library(flowCore)
library(robustbase)
library(ggplot2)
library(ggfortify)
library(GGally)
library(scales)
library(gplots)
library(caret)
library(reshape2)
library(pROC)
library(glmnet)
library(randomForest)
library(glinternet)
library(preprocessCore)
library(FLOWMAPR)
library(grid)

SaveConfusionMatrix <- function(model, data, response, output.name, cutoff) {
  prediction <- predict(model, newdata = data,
                        s = model$bestTune$lambda, type = "prob")
  prediction <- prediction[, 2]
  fitted.prediction <- ifelse(prediction > cutoff, TRUE, FALSE)
  capture.output(confusionMatrix(data = fitted.prediction, data[, response]),
                 file = output.name)
  return()
}

SaveConfusionMatrixGLM <- function(model, data, response, output.name, cutoff) {
  prediction <- predict(model, newdata = data, type = "response")
  fitted.prediction <- ifelse(prediction > cutoff, TRUE, FALSE)
  capture.output(confusionMatrix(data = fitted.prediction, data[, response]),
                 file = output.name)
  return()
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
  svg(paste(base.name, ".svg", sep = ""), width = 4,
      height = 4)
  plot(this.roc, print.auc = TRUE, col = "red", print.auc.y = 5,
       print.auc.x = 85)
  dev.off()
  return(this.roc)
}

SaveROCGLM <- function(model, data, response, output.name,
                       ci, smooth) {
  prediction <- predict(model, newdata = data, type = "response")
  this.roc <- roc(data[, response], prediction, percent = TRUE, algorithm = 1,
                  quiet = TRUE, smooth = smooth, auc = TRUE, ci = ci, plot = TRUE)
  saveRDS(this.roc, file = output.name)
  base.name <- unlist(strsplit(output.name, split = "\\."))[1]
  svg(paste(base.name, ".svg", sep = ""), width = 4,
      height = 4)
  plot(this.roc, print.auc = TRUE, col = "red", print.auc.y = 5,
       print.auc.x = 85)
  dev.off()
  return(this.roc)
}

PlotMultipleROC <- function(rocs, labels, output.name) {
  my.palette <- colorRampPalette(c("red", "orange", "yellow", "green", "cyan", "blue"))
  col <- my.palette(length(rocs))
  svg(paste(output.name, "_auc.svg", sep = ""), width = 4, height = 4)
  for (i in 1:length(rocs)) {
    print(rocs[[i]])
    plot(rocs[[i]], print.auc = TRUE, col = col[i], print.auc.y = 55 - (10 * i),
         add = ifelse(i == 1, FALSE, TRUE))
  }
  dev.off()
  
  svg(paste(output.name, "_legend.svg", sep = ""), width = 4, height = 4)
  for (i in 1:length(rocs)) {
    print(rocs[[i]])
    plot(rocs[[i]], col = col[i], add = ifelse(i == 1, FALSE, TRUE))
  }
  legend("bottomright", legend = labels, col = col, lwd = 5, bty = "n")
  dev.off()
  return()  
}

PlotMultipleROCNoText <- function(rocs, output.name) {
  my.palette <- colorRampPalette(c("red", "orange", "yellow", "green", "cyan", "blue"))
  col <- my.palette(length(rocs))
  svg(paste(output.name, "_auc_no_text.svg", sep = ""), width = 4, height = 4)
  for (i in 1:length(rocs)) {
    print(rocs[[i]])
    plot(rocs[[i]], print.auc = FALSE, col = col[i], xlab = "", ylab = "",
         add = ifelse(i == 1, FALSE, TRUE), labels = FALSE)
  }
  dev.off()
  return()  
}

PlotRFVarImp <- function(rf.model) {
  var.imp <- rf.model$finalModel$importance
  var.imp <- as.data.frame(var.imp[order(var.imp), ])
  barplot(t(var.imp), horiz = TRUE, col = "blue", names = rownames(var.imp), las = 1,
          xlab = "Mean Decrease Gini")
}

PlotGBMVarImp <- function(gbm.model) {
  var.imp <- summary(gbm.model)
  temp <- var.imp
  var.imp <- as.data.frame(var.imp[, 2])
  var.imp <- as.data.frame(var.imp[order(var.imp), ])
  rownames(var.imp) <- rownames(temp)[order(temp[, 2])]
  barplot(t(var.imp), horiz = TRUE, col = "blue", names = rownames(var.imp), las = 1,
          xlab = "Relative Influence")
}

PlotTreeVarImp <- function(model.fit, option = c("gbm", "rf"),
                           output.name) {
  svg(output.name, width = 6, height = 4)
  if (option == "gbm") {
    print(PlotGBMVarImp(model.fit))
  } else if (option == "rf") {
    print(PlotRFVarImp(model.fit))
  } else {
    stop("Do not recognize option!")
  }
  dev.off()
}

PlotTreeVarImpNoText <- function(model.fit, option = c("gbm", "rf"),
                                 output.name) {
  svg(output.name, width = 6, height = 4)
  if (option == "gbm") {
    print(PlotGBMVarImp(model.fit))
  } else if (option == "rf") {
    print(PlotRFVarImpNoText(model.fit))
  } else {
    stop("Do not recognize option!")
  }
  dev.off()
}

PlotRFVarImpNoText <- function(rf.model) {
  var.imp <- rf.model$finalModel$importance
  var.imp <- as.data.frame(var.imp[order(var.imp), ])
  barplot(t(var.imp), horiz = TRUE, col = "blue", names = c(), las = 1,
          xlab = "", names.arg = NULL, col.axis = "white")
}

MakeHeatmap <- function(vals, dendrogram = "column", which.scale = "column",
                        colorpalette, ylab = "Sample", xlab = "Markers") {
  if (!is.null(vals[, ylab])) {
    lab.row <- vals[, ylab]
    vals <- vals[, setdiff(colnames(vals), ylab)]
  } else {
    lab.row <- rownames(vals)
  }
  this.heatmap <- heatmap.2(as.matrix(vals), col = colorpalette,
                            dendrogram = dendrogram,
                            scale = which.scale, ylab = ylab, xlab = xlab,
                            symm = FALSE, density.info = "none", trace = "none",
                            labRow = lab.row, margins = c(6, 6), cexRow = 1, cexCol = 1)
}

MakeHeatmapNoRowReorder <- function(vals, which.scale = "column",
                                    colorpalette, ylab = "Sample", xlab = "Markers") {
  if (!is.null(vals[, ylab])) {
    lab.row <- vals[, ylab]
    vals <- vals[, setdiff(colnames(vals), ylab)]
  } else {
    lab.row <- rownames(vals)
  }
  this.heatmap <- heatmap.2(as.matrix(vals), col = colorpalette,
                            dendrogram = "column", Rowv = FALSE,
                            scale = which.scale, ylab = ylab, xlab = xlab,
                            symm = FALSE, density.info = "none", trace = "none",
                            labRow = lab.row, margins = c(6, 6), cexRow = 1, cexCol = 1)
}

MakeHeatmapNoReorder <- function(vals, which.scale = "column",
                                 colorpalette, ylab = "Sample", xlab = "Markers") {
  if (!is.null(vals[, ylab])) {
    lab.row <- vals[, ylab]
    vals <- vals[, setdiff(colnames(vals), ylab)]
  } else {
    lab.row <- rownames(vals)
  }
  this.heatmap <- heatmap.2(as.matrix(vals), col = colorpalette,
                            dendrogram = "none", Rowv = FALSE, Colv = FALSE,
                            scale = which.scale, ylab = ylab, xlab = xlab,
                            symm = FALSE, density.info = "none", trace = "none",
                            labRow = lab.row, margins = c(6, 6), cexRow = 1, cexCol = 1)
}

MakeHeatmapFlipped <- function(vals, which.scale = "row",
                               colorpalette, xlab = "Sample", ylab = "Markers") {
  if (!is.null(vals[xlab, ])) {
    lab.col <- vals[xlab, ]
    vals <- vals[setdiff(rownames(vals), xlab), ]
  } else {
    lab.col <- colnames(vals)
  }
  save.names <- rownames(vals)
  vals <- as.matrix(sapply(vals, as.numeric))
  rownames(vals) <- save.names
  this.heatmap <- heatmap.2(as.matrix(vals), col = colorpalette,
                            dendrogram = "none", Rowv = FALSE, Colv = FALSE,
                            scale = which.scale, ylab = ylab, xlab = xlab,
                            symm = FALSE, density.info = "none", trace = "none",
                            labCol = lab.col, margins = c(6, 6), cexRow = 1, cexCol = 1)
}

MakeHeatmapFlippedNoText <- function(vals, which.scale = "row",
                               colorpalette, xlab = "Sample", ylab = "Markers") {
  if (!is.null(vals[xlab, ])) {
    lab.col <- vals[xlab, ]
    vals <- vals[setdiff(rownames(vals), xlab), ]
  } else {
    lab.col <- colnames(vals)
  }
  save.names <- rownames(vals)
  vals <- as.matrix(sapply(vals, as.numeric))
  rownames(vals) <- save.names
  this.heatmap <- heatmap.2(as.matrix(vals), col = colorpalette,
                            dendrogram = "none", Rowv = FALSE, Colv = FALSE,
                            scale = which.scale, ylab = "", xlab = "",
                            symm = FALSE, density.info = "none", trace = "none",
                            labRow = FALSE, labCol = FALSE, margins = c(6, 6))
}

MakeHeatmapAsymm <- function(vals, which.scale = "column",
                             colorpalette, ylab = "Sample", xlab = "Markers") {
  if (!is.null(vals[, ylab])) {
    lab.row <- vals[, ylab]
    vals <- vals[, setdiff(colnames(vals), ylab)]
  } else {
    lab.row <- rownames(vals)
  }
  this.heatmap <- heatmap.2(as.matrix(vals), col = colorpalette, symbreaks = FALSE,
                            dendrogram = "none", Rowv = FALSE, Colv = FALSE,
                            scale = which.scale, ylab = ylab, xlab = xlab,
                            symm = FALSE, density.info = "none", trace = "none",
                            symkey = FALSE, breaks = seq(0, 3, length.out = 101),
                            labRow = lab.row, margins = c(6, 6), cexRow = 1, cexCol = 1)
}

MakeHeatmapAsymmNoText <- function(vals, which.scale = "column",
                             colorpalette, ylab = "Sample", xlab = "Markers") {
  if (!is.null(vals[, ylab])) {
    lab.row <- vals[, ylab]
    vals <- vals[, setdiff(colnames(vals), ylab)]
  } else {
    lab.row <- rownames(vals)
  }
  this.heatmap <- heatmap.2(as.matrix(vals), col = colorpalette, symbreaks = FALSE,
                            dendrogram = "none", Rowv = FALSE, Colv = FALSE,
                            scale = which.scale, ylab = "", xlab = "",
                            symm = FALSE, density.info = "none", trace = "none",
                            symkey = FALSE, breaks = seq(0, 3, length.out = 101),
                            labRow = FALSE, labCol = FALSE, margins = c(6, 6))
}

rm(list = ls())

pre.folder <- "/Users/mesako/Desktop/20180211_Bcl-2_Analysis/"
setwd(pre.folder)
source("20180710_Bcl-2_Figures_Support.R")

seed.X <- 1
set.seed(seed.X)
analysis.mode <- "lasttime"
pre <- analysis.mode

############################
### Establish Parameters ###
############################

signaling.channels <- c("Akt_pS473", "p53", "STAT5_pY694",
                        "cMyc", "IkBa", "S6_pS235", "CREB_pS133",
                        "p38_pT180", "Erk_p44")
pro.survival.Bcl2 <- c("Mcl1", "Bcl2", "BclxL")
pro.apoptotic.Bcl2 <- c("Bim", "Bax", "Bak_activated")
Bcl2.family.channels <- c(pro.survival.Bcl2, pro.apoptotic.Bcl2)
most.channels <- c(Bcl2.family.channels, signaling.channels)

setwd(paste(pre.folder, "run3_live_classification_models", sep = ""))
setwd(paste("run3_live_", pre,"_classification_models", sep = ""))

##################################
### Reproduce Modeling Figures ###
##################################

Bcl2.rf.bor.fit <- readRDS(file = paste("models/", pre, "_live_Bcl2_bor_rf_fit_run3.rds", sep = ""))
Bcl2.rf.dex.fit <- readRDS(file = paste("models/", pre, "_live_Bcl2_dex_rf_fit_run3.rds", sep = ""))
most.rf.bor.fit <- readRDS(file = paste("models/", pre, "_live_most_bor_rf_fit_run3.rds", sep = ""))
most.rf.dex.fit <- readRDS(file = paste("models/", pre, "_live_most_dex_rf_fit_run3.rds", sep = ""))

most.rf.bor.roc <- readRDS(paste("roc_curves/", pre, "_live_most_rf_bor_roc.rds", sep = ""))
most.rf.dex.roc <- readRDS(paste("roc_curves/", pre, "_live_most_rf_dex_roc.rds", sep = ""))
Bcl2.rf.bor.roc <- readRDS(paste("roc_curves/", pre, "_live_Bcl2_rf_bor_roc.rds", sep = ""))
Bcl2.rf.dex.roc <- readRDS(paste("roc_curves/", pre, "_live_Bcl2_rf_dex_roc.rds", sep = ""))

most.rf.bor.on.dex.roc <- readRDS(paste("roc_curves/", pre, "_live_most_rf_bor_on_dex_roc.rds", sep = ""))
most.rf.dex.on.bor.roc <- readRDS(paste("roc_curves/", pre, "_live_most_rf_dex_on_bor_roc.rds", sep = ""))
Bcl2.rf.bor.on.dex.roc <- readRDS(paste("roc_curves/", pre, "_live_Bcl2_rf_bor_on_dex_roc.rds", sep = ""))
Bcl2.rf.dex.on.bor.roc <- readRDS(paste("roc_curves/", pre, "_live_Bcl2_rf_dex_on_bor_roc.rds", sep = ""))

most.rf.bor.on.DMSO.roc <- readRDS(paste("roc_curves/", pre, "_live_most_rf_bor_on_DMSO_roc.rds", sep = ""))
most.rf.dex.on.DMSO.roc <- readRDS(paste("roc_curves/", pre, "_live_most_rf_dex_on_DMSO_roc.rds", sep = ""))
Bcl2.rf.bor.on.DMSO.roc <- readRDS(paste("roc_curves/", pre, "_live_Bcl2_rf_bor_on_DMSO_roc.rds", sep = ""))
Bcl2.rf.dex.on.DMSO.roc <- readRDS(paste("roc_curves/", pre, "_live_Bcl2_rf_dex_on_DMSO_roc.rds", sep = ""))

run1.rf.most.bor.roc <- readRDS(paste("roc_curves/", pre, "_live_run1_rf_most_bor_roc.rds", sep = ""))
run2.rf.most.bor.roc <- readRDS(paste("roc_curves/", pre, "_live_run2_rf_most_bor_roc.rds", sep = ""))
run1.rf.most.dex.roc <- readRDS(paste("roc_curves/", pre, "_live_run1_rf_most_dex_roc.rds", sep = ""))
run2.rf.most.dex.roc <- readRDS(paste("roc_curves/", pre, "_live_run2_rf_most_dex_roc.rds", sep = ""))

run1.rf.Bcl2.bor.roc <- readRDS(paste("roc_curves/", pre, "_live_run1_rf_Bcl2_bor_roc.rds", sep = ""))
run2.rf.Bcl2.bor.roc <- readRDS(paste("roc_curves/", pre, "_live_run2_rf_Bcl2_bor_roc.rds", sep = ""))
run1.rf.Bcl2.dex.roc <- readRDS(paste("roc_curves/", pre, "_live_run1_rf_Bcl2_dex_roc.rds", sep = ""))
run2.rf.Bcl2.dex.roc <- readRDS(paste("roc_curves/", pre, "_live_run2_rf_Bcl2_dex_roc.rds", sep = ""))

most.rf.dex.on.BakBaxDKO.roc <- readRDS(paste("roc_curves/", pre, "_live_accuracy_most_rf_dex_fit_on_BakBaxDKO_roc.rds", sep = ""))
Bcl2.rf.dex.on.BakBaxDKO.roc <- readRDS(paste("roc_curves/", pre, "_live_accuracy_Bcl2_rf_dex_fit_on_BakBaxDKO_roc.rds", sep = ""))

##########################
### Make FINAL Figures ###
##########################

setwd("/Users/mesako/Desktop/20180211_Bcl-2_Analysis/20180710_Bcl-2_New_SVG_Figures")

# FIGURE 1C
save.palette <- colorRampPalette(c("#1500FB", "#C3C3C7", "#D10100"))
save.palette <- save.palette(100)

asymm.palette1 <- colorRampPalette(c("#1500FB", "#C3C3C7"))
asymm.palette2 <- colorRampPalette(c("#C3C3C7", "#D10100"))
asymm.palette <- asymm.palette1(34)
asymm.palette <- c(asymm.palette, asymm.palette2(67)[-1])

cell.line.data <- read.csv("/Users/mesako/Downloads/Bcl2_WB_CyTOF_Comparison/20180430_MM_cell_lines_Bcl2_mean_stats.csv", header = TRUE)

colnames(cell.line.data) <- c("FCS_File", "Mcl1", "Bcl2", "BclxL", "Bim",
                              "active_Bak", "Bak", "active_Bax", "Bax")
cell.line.data$FCS_File <- c("AMO1", "EJM", "H929", "KMS-12-BM", "KMS-12-PE",
                             "KMS-28-BM", "KMS-28-PE", "KMS-34", "MM1R",
                             "MM1S", "OPM2", "U266B1")
cell.line.data[2:ncol(cell.line.data)] <- apply(cell.line.data[2:ncol(cell.line.data)], 2, FLOWMAPR:::Asinh)
cell.line.data <- cell.line.data[, setdiff(colnames(cell.line.data), c("active_Bax", "Bak"))]
cell.line.data <- cell.line.data[, c(1:5, 7, 6)]

first.row <- t(cell.line.data)[1, ]
flip.data <- as.data.frame(t(cell.line.data[, 2:ncol(cell.line.data)]))
flip.data <- rbind(FCS_File = first.row, flip.data)
for (x in 2:nrow(flip.data)) {
  flip.data[x, ] <- as.numeric(flip.data[x, ])
}

plot.name <- paste("MM_cell_lines_Bcl2_profiles_flipped_heatmap_no_text.svg", sep = "")
# plot.name <- paste("MM_cell_lines_Bcl2_profiles_flipped_heatmap.svg", sep = "")
svg(plot.name, width = 8, height = 6)
MakeHeatmapFlippedNoText(flip.data, which.scale = "row",
                         colorpalette = save.palette, xlab = "FCS_File", ylab = "Markers")
# MakeHeatmapFlipped(flip.data, which.scale = "row",
#                    colorpalette = save.palette, xlab = "FCS_File", ylab = "Markers")
dev.off()

# FIGURE 1E
WB.data <- read.csv("/Users/mesako/Downloads/Bcl2_WB_CyTOF_Comparison/20180213_MM_cell_lines_WB_data.csv", header = TRUE)
CyTOF.data <- read.csv("/Users/mesako/Downloads/Bcl2_WB_CyTOF_Comparison/20180430_MM_cell_lines_Bcl2_mean_stats.csv", header = TRUE)
WB.data$Cell_Line <- as.character(WB.data$Cell_Line)
CyTOF.data$FCS.Filename <- as.character(CyTOF.data$FCS.Filename)

WB.data <- WB.data[!(WB.data$Cell_Line == ""), ]
WB.data <- WB.data[, setdiff(colnames(WB.data), "HSP70")]
colnames(WB.data)[grepl(colnames(WB.data), pattern = "Bak")] <- "active_Bak"
WB.data <- WB.data[, c(1:4, 7, 5:6)]

colnames(CyTOF.data) <- c("Cell_Line", "active_Bax", "Bak", "BclxL", "Bax",
                          "active_Bak", "Bcl2", "Mcl1", "Bim")
CyTOF.data$Cell_Line <- c("AMO1", "EJM", "H929", "KMS-12-BM", "KMS-12-PE",
                          "KMS-28-BM", "KMS-28-PE", "KMS-34", "MM1R",
                          "MM1S", "OPM2", "U266B1")

CyTOF.data[2:ncol(CyTOF.data)] <- apply(CyTOF.data[2:ncol(CyTOF.data)], 2, FLOWMAPR:::Asinh)
CyTOF.data <- CyTOF.data[, setdiff(colnames(CyTOF.data), c("active_Bax", "Bak"))]
CyTOF.data[2:ncol(CyTOF.data)] <- apply(CyTOF.data[2:ncol(CyTOF.data)], 2, rescale)
CyTOF.data <- CyTOF.data[, c(1, 6, 5, 2, 7, 3, 4)]

for (x in 2:ncol(CyTOF.data)) {
  this.compare <- colnames(CyTOF.data)[x]
  v1 <- CyTOF.data[, this.compare]
  v2 <- WB.data[, this.compare]
  temp <- data.frame(v1, v2)
  this.cor1 <- cor(v1, v2, method = "pearson")
  this.cor2 <- cor(v1, v2, method = "spearman")
  this.cor1 <- round(this.cor1, digits = 4)
  this.cor2 <- round(this.cor2, digits = 4)
  cat("CyTOF", this.compare, "and WB", this.compare, "Pearson","=", this.cor1, "\n")
  cat("CyTOF", this.compare, "and WB", this.compare, "Spearman", "=", this.cor2, "\n")
  grob <- grobTree(textGrob(paste("Pearson Cor:", this.cor1, "\n", "Spearman Cor:", this.cor2, sep = " "),
                            x = 0.1,  y = 0.95, hjust = 0,
                            gp = gpar(col="black", fontsize = 10)))
  this.plot <- ggplot(temp, aes(v1, v2)) + geom_point() +  
    geom_smooth(method = "lm") + xlab(paste("CyTOF_", this.compare, sep = "")) + 
    ylab(paste("WB_", this.compare, sep = ""))
  # + geom_text(aes(label = CyTOF.data[, 1]), hjust = 0.5, vjust = -0.5)
  # this.plot <- this.plot + annotation_custom(grob)
  this.plot <- this.plot + theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())
  
  plot.name <- paste("CyTOF_vs_WB_", this.compare, "_no_text.svg", sep = "")
  # plot.name <- paste("CyTOF_vs_WB_", this.compare, ".svg", sep = "")
  svg(plot.name, width = 6, height = 6)
  print(this.plot)
  dev.off()
}

WB.data[2:ncol(WB.data)] <- apply(WB.data[2:ncol(WB.data)], 2, rescale)

for (x in 2:ncol(CyTOF.data)) {
  this.compare <- colnames(CyTOF.data)[x]
  v1 <- CyTOF.data[, this.compare]
  v2 <- WB.data[, this.compare]
  temp <- data.frame(v1, v2)
  this.cor1 <- cor(v1, v2, method = "pearson")
  this.cor2 <- cor(v1, v2, method = "spearman")
  this.cor1 <- round(this.cor1, digits = 4)
  this.cor2 <- round(this.cor2, digits = 4)
  cat("CyTOF", this.compare, "and WB", this.compare, "Pearson","=", this.cor1, "\n")
  cat("CyTOF", this.compare, "and WB", this.compare, "Spearman", "=", this.cor2, "\n")
  grob <- grobTree(textGrob(paste("Pearson Cor:", this.cor1, "\n", "Spearman Cor:", this.cor2, sep = " "),
                            x = 0.1,  y = 0.95, hjust = 0,
                            gp = gpar(col="black", fontsize = 10)))
  this.plot <- ggplot(temp, aes(v1, v2)) + geom_point() +  
    geom_smooth(method = "lm") + xlab(paste("CyTOF_", this.compare, sep = "")) + 
    ylab(paste("WB_", this.compare, sep = ""))
  # + geom_text(aes(label = CyTOF.data[, 1]), hjust = 0.5, vjust = -0.5)
  # this.plot <- this.plot + annotation_custom(grob)
  this.plot <- this.plot + theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())  
  plot.name <- paste("CyTOF_vs_WB_", this.compare, "_WB_rescaled_no_text.svg", sep = "")
  # plot.name <- paste("CyTOF_vs_WB_", this.compare, "_WB_rescaled.svg", sep = "")
  svg(plot.name, width = 6, height = 6)
  print(this.plot)
  dev.off()
}

# FIGURE 2B
treatment.data1 <- read.csv("/Users/mesako/Downloads/Bcl2_Heatmaps/20180328_MM1S_run_3A_exported_stats.csv", header = TRUE)
treatment.data2 <- read.csv("/Users/mesako/Downloads/Bcl2_Heatmaps/20180328_MM1S_run_3B_exported_stats.csv", header = TRUE)
treatment.data <- rbind(treatment.data1, treatment.data2)

treatment.data <- treatment.data[!grepl(treatment.data$FCS.Filename, pattern = "HELA"), ]
treatment.data <- treatment.data[grepl(treatment.data$FCS.Filename, pattern = "WT"), ]
treatment.data <- treatment.data[!grepl(treatment.data$FCS.Filename, pattern = "DMSO"), ]
treatment.data$FCS.Filename <- as.character(treatment.data$FCS.Filename)
treatment.data$FCS.Filename <- c("Untreated", "iBclxL_01", "iBclxL_06", "Len_06",
                                 "Len_24", "Len_48", "Len_72", "Dex_24", "Dex_72",
                                 "JQ1_24", "Bor_06", "Bor_24")

treatment.data[2:ncol(treatment.data)] <- apply(treatment.data[2:ncol(treatment.data)], 2, FLOWMAPR:::Asinh)
save.treat.data <- treatment.data

singlets.col <- colnames(treatment.data)[grepl(colnames(treatment.data), pattern = "Singlets")]

fix.col <- c("p-p38", "pBcl2", "pErk", "pRb", "pAkt", "activeBak",
             "CyclinB", "pSTAT5", "CyclinA", "pH3", "pBadS112",
             "pZAP70", "pS6", "pCREB")
remove.markers <- c("pBcl2", "APAF", "Bclw", "pZAP70", "pBadS112", "Bak")
cell.cycle.markers <- c("CyclinB", "pRb", "CyclinA", "pH3")
marker.order <- c("Sample", "Cisplatin", "aCasp3", "cPARP", "Bcl2",
                  "BclxL", "Mcl1", "Bim", "Bax", "activeBak",
                  "pH3", "CyclinA", "CyclinB", "pRb", "pCREB",
                  "p-p38", "IkBa","pSTAT5", "pS6", "pAkt", "pErk",
                  "p53", "cMyc")

singlets.ratio.data <- save.treat.data[, c("FCS.Filename", singlets.col)]
temp <- unlist(strsplit(colnames(singlets.ratio.data), split = "\\.."))
temp <- temp[seq(from = 6, to = length(temp), by = 8)]
temp <- c("Sample", temp)
colnames(singlets.ratio.data) <- temp
colnames(singlets.ratio.data)[grepl(colnames(singlets.ratio.data), pattern = "\\_")] <- fix.col

singlets.ratio.data <- singlets.ratio.data[, setdiff(colnames(singlets.ratio.data), remove.markers)]
singlets.ratio.data <- singlets.ratio.data[, marker.order]
for (x in 2:ncol(singlets.ratio.data)) {
  singlets.ratio.data[, x] <- singlets.ratio.data[, x] / singlets.ratio.data[1, x]
}

temp.ratio.data <- singlets.ratio.data[c(1, 11, 12, 8, 9), ]
plot.name <- paste("MM1S_treatment_bor_dex_only_singlets_profiles_asinh_medians_ratio_heatmap_no_reorder_no_text.svg", sep = "")
# plot.name <- paste("MM1S_treatment_bor_dex_only_singlets_profiles_asinh_medians_ratio_heatmap_no_reorder.svg", sep = "")
svg(plot.name, width = 10, height = 5)
MakeHeatmapAsymmNoText(temp.ratio.data, which.scale = "none",
                       colorpalette = asymm.palette, ylab = "Sample", xlab = "Markers")
# MakeHeatmapAsymm(temp.ratio.data, which.scale = "none",
#                  colorpalette = asymm.palette, ylab = "Sample", xlab = "Markers")
dev.off()

# FIGURE 4B
# PlotTreeVarImp(most.rf.bor.fit, "rf", paste(pre, "_most_rf_bor_fit_varimp_run3.svg", sep = ""))
# PlotTreeVarImp(most.rf.dex.fit, "rf", paste(pre, "_most_rf_dex_fit_varimp_run3.svg", sep = ""))

PlotTreeVarImpNoText(most.rf.bor.fit, "rf", paste(pre, "_most_rf_bor_fit_varimp_run3_no_text.svg", sep = ""))
PlotTreeVarImpNoText(most.rf.dex.fit, "rf", paste(pre, "_most_rf_dex_fit_varimp_run3_no_text.svg", sep = ""))

# FIGURE 4C
# PlotTreeVarImp(Bcl2.rf.bor.fit, "rf", paste(pre, "_Bcl2_rf_bor_fit_varimp_run3.svg", sep = ""))
# PlotTreeVarImp(Bcl2.rf.dex.fit, "rf", paste(pre, "_Bcl2_rf_dex_fit_varimp_run3.svg", sep = ""))

PlotTreeVarImpNoText(Bcl2.rf.bor.fit, "rf", paste(pre, "_Bcl2_rf_bor_fit_varimp_run3_no_text.svg", sep = ""))
PlotTreeVarImpNoText(Bcl2.rf.dex.fit, "rf", paste(pre, "_Bcl2_rf_dex_fit_varimp_run3_no_text.svg", sep = ""))

# FIGURE 5B
# PlotMultipleROC(list(run1.rf.most.bor.roc, run2.rf.most.bor.roc, most.rf.bor.roc),
#                 c("Run1", "Run2", "Run3"),
#                 output.name = "run_comparison_rf_bortezomib_most_models")
# PlotMultipleROC(list(run1.rf.most.dex.roc, run2.rf.most.dex.roc, most.rf.dex.roc),
#                 c("Run1", "Run2", "Run3"),
#                 output.name = "run_comparison_rf_dexamethasone_most_models")
# 
# PlotMultipleROC(list(run1.rf.Bcl2.bor.roc, run2.rf.Bcl2.bor.roc, Bcl2.rf.bor.roc),
#                 c("Run1", "Run2", "Run3"),
#                 output.name = "run_comparison_rf_bortezomib_Bcl2_models")
# PlotMultipleROC(list(run1.rf.Bcl2.dex.roc, run2.rf.Bcl2.dex.roc, Bcl2.rf.dex.roc),
#                 c("Run1", "Run2", "Run3"),
#                 output.name = "run_comparison_rf_dexamethasone_Bcl2_models")

PlotMultipleROCNoText(list(run1.rf.most.bor.roc, run2.rf.most.bor.roc, most.rf.bor.roc),
                      output.name = "run_comparison_rf_bortezomib_most_models")
PlotMultipleROCNoText(list(run1.rf.most.dex.roc, run2.rf.most.dex.roc, most.rf.dex.roc),
                      output.name = "run_comparison_rf_dexamethasone_most_models")

PlotMultipleROCNoText(list(run1.rf.Bcl2.bor.roc, run2.rf.Bcl2.bor.roc, Bcl2.rf.bor.roc),
                      output.name = "run_comparison_rf_bortezomib_Bcl2_models")
PlotMultipleROCNoText(list(run1.rf.Bcl2.dex.roc, run2.rf.Bcl2.dex.roc, Bcl2.rf.dex.roc),
                      output.name = "run_comparison_rf_dexamethasone_Bcl2_models")

# FIGURE 6D
# PlotMultipleROC(list(most.rf.dex.roc, most.rf.dex.on.BakBaxDKO.roc),
#                 c("WT", "BakBaxDKO"),
#                 output.name = "live_most_rf_dexamethasone_on_BakBaxDKO")
# 
# PlotMultipleROC(list(Bcl2.rf.dex.roc, Bcl2.rf.dex.on.BakBaxDKO.roc),
#                 c("WT", "BakBaxDKO"),
#                 output.name = "live_Bcl2_rf_dexamethasone_on_BakBaxDKO")

PlotMultipleROCNoText(list(most.rf.dex.roc, most.rf.dex.on.BakBaxDKO.roc),
                      output.name = "live_most_rf_dexamethasone_on_BakBaxDKO")

PlotMultipleROCNoText(list(Bcl2.rf.dex.roc, Bcl2.rf.dex.on.BakBaxDKO.roc),
                      output.name = "live_Bcl2_rf_dexamethasone_on_BakBaxDKO")

# FIGURE 7C
# PlotMultipleROC(list(most.rf.bor.roc, most.rf.dex.roc,
#                      most.rf.bor.on.dex.roc,
#                      most.rf.dex.on.bor.roc),
#                 c("Bor on Bor", "Dex on Dex", "Bor on Dex", "Dex on Bor"),
#                 output.name = "live_most_rf_cross-drug_comparison_models")
# 
# PlotMultipleROC(list(Bcl2.rf.bor.roc, Bcl2.rf.dex.roc,
#                      Bcl2.rf.bor.on.dex.roc,
#                      Bcl2.rf.dex.on.bor.roc),
#                 c("Bor on Bor", "Dex on Dex", "Bor on Dex", "Dex on Bor"),
#                 output.name = "live_Bcl2_rf_cross-drug_comparison_models")

PlotMultipleROCNoText(list(most.rf.bor.roc, most.rf.dex.roc,
                           most.rf.bor.on.dex.roc,
                           most.rf.dex.on.bor.roc),
                      output.name = "live_most_rf_cross-drug_comparison_models")

PlotMultipleROCNoText(list(Bcl2.rf.bor.roc, Bcl2.rf.dex.roc,
                           Bcl2.rf.bor.on.dex.roc,
                           Bcl2.rf.dex.on.bor.roc),
                      output.name = "live_Bcl2_rf_cross-drug_comparison_models")

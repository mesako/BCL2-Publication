rm(list = ls())

pre.folder <- "/Users/mesako/Desktop/20180211_Bcl-2_Analysis/"
setwd(pre.folder)
source("figures_support.R")

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

# FIGURE 4E
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

# FIGURE 5E
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

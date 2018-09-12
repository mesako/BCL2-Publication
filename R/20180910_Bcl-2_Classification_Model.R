rm(list = ls())

pre.folder <- "/Users/mesako/Desktop/20180211_Bcl-2_Analysis/"
setwd(pre.folder)
source("20180910_Bcl-2_Classification_Model_Support.R")

seed.X <- 1
set.seed(seed.X)

# analysis.mode <- "alltime"
analysis.mode <- "lasttime"
pre <- analysis.mode

###################################################
### Load All FCS Files and Establish Parameters ###
###################################################

all.files.folder <- "/Users/mesako/Desktop/20180211_Bcl-2_Analysis/FCS_Files_Live_Together"
all.file.names <- FLOWMAPR::GetFCSNames(all.files.folder)
var.annotate <- FLOWMAPR::ConstructVarAnnotate(all.file.names[1])
var.remove <- FLOWMAPR::SuggestVarRemove(var.annotate)
var.remove <- c(var.remove, "File_Number", "Yb171Di", "CD61",
                "ZAP70_pY319", "Bi209Di")
file.names.short <- basename(all.file.names)
file.names.short <- gsub(file.names.short, pattern = "\\.fcs", replacement = "")

###################################
### Quantile Normalize All Data ###
###################################

all.file.names <- all.file.names[1]
all.files.together <- FLOWMAPR::LoadCleanFCS(all.file.names, var.remove,
                                             var.annotate, subsamples = FALSE)
names(all.files.together) <- file.names.short
saveRDS(all.files.together, file = "all_runs_live_cleaned.rds")

### checking aCasp3 levels before and after normalization
### in dexamethasone samples, last timepoint (72 hours)
all.files.together <- readRDS(file = "all_runs_live_cleaned.rds")

all.files.together.norm <- all.files.together
for (i in 1:length(all.files.together.norm)) {
  temp <- as.matrix(all.files.together.norm[[i]])
  temp <- normalize.quantiles(temp)
  temp <- as.data.frame(temp)
  colnames(temp) <- colnames(all.files.together[[i]])
  all.files.together.norm[[i]] <- temp
}
rm(all.files.together)

names(all.files.together.norm) <- file.names.short
saveRDS(all.files.together.norm, file = "all_runs_live_cleaned_quant_norm.rds")
rm(all.files.together.norm)

#######################################################
### Load Quantile Normalized or Non-Normalized Data ###
#######################################################

# Use Run 3, check aCasp3 bimodality on singlets gated cells
setwd(pre.folder)
which.use <- "no_norm"
# which.use <- "norm"
if (which.use == "norm") {
  print("use normed")
  all.files <- readRDS("all_runs_live_cleaned_quant_norm.rds")
} else if (which.use == "no_norm") {
  print("use non-normed")
  all.files <- readRDS("all_runs_live_cleaned.rds")
} else {
  stop("not recognized")
}
files.to.use <- names(all.files)[grepl(names(all.files), pattern = "RunThree")]
files.to.use <- files.to.use[grepl(files.to.use, pattern = "WT")]

bor.files <- c(files.to.use[grepl(files.to.use, pattern = "untreated")],
               files.to.use[grepl(files.to.use, pattern = "B-")])
dex.files <- c(files.to.use[grepl(files.to.use, pattern = "untreated")],
               files.to.use[grepl(files.to.use, pattern = "D-")])
bor.data <- all.files[bor.files]
dex.data <- all.files[dex.files]
bor.data <- CollapseBy(bor.data, c("00", "06", "24"),
                       "Timepoint", subsample = FALSE)
dex.data <- CollapseBy(dex.data, c("00", "24", "72"),
                       "Timepoint", subsample = FALSE)
rm(all.files)

##########################################################
### Plot aCasp3 Boundaries for Apoptotic/Non-Apoptotic ###
##########################################################

setwd(paste(pre.folder, "run3_live_classification_models", sep = ""))

temp.bor.data <- bor.data[bor.data[, "Timepoint"] == "24", ]
temp.dex.data <- dex.data[dex.data[, "Timepoint"] == "72", ]

if (which.use == "norm") {
  # Normalized Values
  bor.low.cutoff <- 2.4
  bor.high.cutoff <- 2.7
  dex.low.cutoff <- 1.3
  dex.high.cutoff <- 1.6
  xlim <- c(0, 6)
} else if (which.use == "no_norm") {
  # Non-Normalized Values
  bor.low.cutoff <- 2.25
  bor.high.cutoff <- 3.25
  dex.low.cutoff <- 2
  dex.high.cutoff <- 3
  xlim <- c(-1, 7)
} else {
  stop("not recognized")
}

qplot(temp.bor.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Bortezomib",
      xlim = xlim) +
  geom_vline(xintercept = bor.low.cutoff, col = "red") + 
  geom_vline(xintercept = bor.high.cutoff, col = "red") 
output.name <- paste(which.use, "_bortezomib_live_lasttime_cCaspase3.png", sep = "")
ggsave(output.name, width = 6, height = 4, dpi = 400)

qplot(bor.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Bortezomib",
      xlim = xlim) +
  geom_vline(xintercept = bor.low.cutoff, col = "red") + 
  geom_vline(xintercept = bor.high.cutoff, col = "red") 
output.name <- paste(which.use, "_bortezomib_live_alltime_cCaspase3.png", sep = "")
ggsave(output.name, width = 6, height = 4, dpi = 400)

qplot(temp.dex.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Dexamethasone",
      xlim = xlim) +
  geom_vline(xintercept = dex.low.cutoff, col = "red") + 
  geom_vline(xintercept = dex.high.cutoff, col = "red") 
output.name <- paste(which.use, "_dexamethasone_live_lasttime_cCaspase3.png", sep = "")
ggsave(output.name, width = 6, height = 4, dpi = 400)

qplot(dex.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Dexamethasone",
      xlim = xlim) +
  geom_vline(xintercept = dex.low.cutoff, col = "red") + 
  geom_vline(xintercept = dex.high.cutoff, col = "red") 
output.name <- paste(which.use, "_dexamethasone_live_alltime_cCaspase3.png", sep = "")
ggsave(output.name, width = 6, height = 4, dpi = 400)

######################################################################
### Gate Out Apoptotic/Non-Apoptotic Cells Using Non-Normed Values ###
######################################################################

setwd(pre.folder)
all.files <- readRDS("all_runs_live_cleaned_quant_norm.rds")
norm.files.to.use <- names(all.files)[grepl(names(all.files), pattern = "RunThree")]
norm.files.to.use <- norm.files.to.use[grepl(norm.files.to.use, pattern = "WT")]
norm.bor.files <- c(norm.files.to.use[grepl(norm.files.to.use, pattern = "untreated")],
                    norm.files.to.use[grepl(norm.files.to.use, pattern = "B-")])
norm.dex.files <- c(norm.files.to.use[grepl(norm.files.to.use, pattern = "untreated")],
                    norm.files.to.use[grepl(norm.files.to.use, pattern = "D-")])
norm.bor.data <- all.files[norm.bor.files]
norm.dex.data <- all.files[norm.dex.files]
rm(all.files)
norm.bor.data <- CollapseBy(norm.bor.data, c("00", "06", "24"),
                            "Timepoint", subsample = FALSE)
norm.dex.data <- CollapseBy(norm.dex.data, c("00", "24", "72"),
                            "Timepoint", subsample = FALSE)

all.files <- readRDS("all_runs_live_cleaned.rds")
files.to.use <- names(all.files)[grepl(names(all.files), pattern = "RunThree")]
files.to.use <- files.to.use[grepl(files.to.use, pattern = "WT")]

bor.files <- c(files.to.use[grepl(files.to.use, pattern = "untreated")],
               files.to.use[grepl(files.to.use, pattern = "B-")])
dex.files <- c(files.to.use[grepl(files.to.use, pattern = "untreated")],
               files.to.use[grepl(files.to.use, pattern = "D-")])
bor.data <- all.files[bor.files]
dex.data <- all.files[dex.files]
rm(all.files)
bor.data <- CollapseBy(bor.data, c("00", "06", "24"),
                       "Timepoint", subsample = FALSE)
dex.data <- CollapseBy(dex.data, c("00", "24", "72"),
                       "Timepoint", subsample = FALSE)

# Non-Normalized Values
bor.low.cutoff <- 2.25
bor.high.cutoff <- 3.25
dex.low.cutoff <- 2
dex.high.cutoff <- 3

which.ind <- which(bor.data$aCasp3 < bor.low.cutoff)
nonapoptotic.bor.data <- norm.bor.data[which.ind, ]
which.ind <- which(dex.data$aCasp3 < dex.low.cutoff)
nonapoptotic.dex.data <- norm.dex.data[which.ind, ]
which.ind <- which(bor.data$aCasp3 > bor.high.cutoff)
apoptotic.bor.data <- norm.bor.data[which.ind, ]
which.ind <- which(dex.data$aCasp3 > dex.high.cutoff)
apoptotic.dex.data <- norm.dex.data[which.ind, ]

#####################################################################
### Recombined Non-Apoptotic/Apoptotic Cells With Status Labeling ###
#####################################################################

setwd("run3_live_classification_models")
nonapoptotic.bor.data <- cbind(nonapoptotic.bor.data, Status = rep("Resistant", times = nrow(nonapoptotic.bor.data)))
nonapoptotic.dex.data <- cbind(nonapoptotic.dex.data, Status = rep("Resistant", times = nrow(nonapoptotic.dex.data)))
apoptotic.bor.data <- cbind(apoptotic.bor.data, Status = rep("Sensitive", times = nrow(apoptotic.bor.data)))
apoptotic.dex.data <- cbind(apoptotic.dex.data, Status = rep("Sensitive", times = nrow(apoptotic.dex.data)))

recombined.bor.data <- rbind(nonapoptotic.bor.data, apoptotic.bor.data)
recombined.dex.data <- rbind(nonapoptotic.dex.data, apoptotic.dex.data)
recombined.bor.data <- recombined.bor.data[sample(1:nrow(recombined.bor.data), size = nrow(recombined.bor.data)), ]
recombined.dex.data <- recombined.dex.data[sample(1:nrow(recombined.dex.data), size = nrow(recombined.dex.data)), ]

recombined.bor.data <- ConvertClassto01(recombined.bor.data, column = "Status", true.label = "Resistant")
recombined.dex.data <- ConvertClassto01(recombined.dex.data, column = "Status", true.label = "Resistant")

saveRDS(recombined.bor.data, file = "run3_bortezomib_data.rds")
saveRDS(recombined.dex.data, file = "run3_dexamethasone_data.rds")

#################################################
### Prepare Data For Modeling By Partitioning ###
#################################################

setwd(paste(pre.folder, "run3_live_classification_models", sep = ""))
recombined.bor.data <- readRDS(file = "run3_bortezomib_data.rds")
recombined.dex.data <- readRDS(file = "run3_dexamethasone_data.rds")

if (analysis.mode == "alltime") {
  bor.data <- recombined.bor.data
  dex.data <- recombined.dex.data
} else if (analysis.mode == "lasttime") {
  bor.data <- recombined.bor.data[recombined.bor.data[, "Timepoint"] == "24", ]
  dex.data <- recombined.dex.data[recombined.dex.data[, "Timepoint"] == "72", ]
} else {
  stop("mode not recognized!")
}

setwd(paste("run3_live_", pre,"_classification_models", sep = ""))

train.ind <- createDataPartition(bor.data$Resistant, p = 0.8, list = FALSE)
training.bor.data <- bor.data[train.ind, ]
test.bor.data <- bor.data[-train.ind, ]
training.bor.data[, "Resistant"] <- as.factor(training.bor.data[, "Resistant"])
test.bor.data[, "Resistant"] <- as.factor(test.bor.data[, "Resistant"])
saveRDS(training.bor.data, file = paste("data/", pre, "_training_live_bortezomib_data_run3.rds", sep = ""))
saveRDS(test.bor.data, file = paste("data/", pre, "_test_live_bortezomib_data_run3.rds", sep = ""))

train.ind <- createDataPartition(dex.data$Resistant, p = 0.8, list = FALSE)
training.dex.data <- dex.data[train.ind, ]
test.dex.data <- dex.data[-train.ind, ]
training.dex.data[, "Resistant"] <- as.factor(training.dex.data[, "Resistant"])
test.dex.data[, "Resistant"] <- as.factor(test.dex.data[, "Resistant"])
saveRDS(training.dex.data, file = paste("data/", pre, "_training_live_dexamethasone_data_run3.rds", sep = ""))
saveRDS(test.dex.data, file = paste("data/", pre, "_test_live_dexamethasone_data_run3.rds", sep = ""))

###########################################
### Basic Modeling For Resistant Status ###
###########################################

setwd(paste("run3_live_", pre,"_classification_models", sep = ""))
training.bor.data <- readRDS(paste("data/", pre, "_training_live_bortezomib_data_run3.rds", sep = ""))
training.dex.data <- readRDS(paste("data/", pre, "_training_live_dexamethasone_data_run3.rds", sep = ""))

signaling.channels <- c("Akt_pS473", "p53", "STAT5_pY694",
                        "cMyc", "IkBa", "S6_pS235", "CREB_pS133",
                        "p38_pT180", "Erk_p44")
pro.survival.Bcl2 <- c("Mcl1", "Bcl2", "BclxL")
pro.apoptotic.Bcl2 <- c("Bim", "Bax", "Bak_activated")
Bcl2.family.channels <- c(pro.survival.Bcl2, pro.apoptotic.Bcl2)
most.channels <- c(Bcl2.family.channels, signaling.channels)

most.bor.data <- subset(training.bor.data, select = c(most.channels, "Resistant"))
most.dex.data <- subset(training.dex.data, select = c(most.channels, "Resistant"))
signaling.bor.data <- subset(training.bor.data, select = c(signaling.channels, "Resistant"))
signaling.dex.data <- subset(training.dex.data, select = c(signaling.channels, "Resistant"))
Bcl2.bor.data <- subset(training.bor.data, select = c(Bcl2.family.channels, "Resistant"))
Bcl2.dex.data <- subset(training.dex.data, select = c(Bcl2.family.channels, "Resistant"))

##################################################################
### Learn Classification Models for Resistant on Training Data ###
######### Use 10-fold Cross-Validation, Repeated 10 Times ########
##################################################################

fit.control <- trainControl(method = "repeatedcv", number = 10, repeats = 10)

most.bor.fit <- train(Resistant ~ .,  data = most.bor.data, method = "glmnet",
                      trControl = fit.control, tuneLength = 5)
saveRDS(most.bor.fit, file = paste("models/", pre, "_live_most_bor_fit_run3.rds", sep = ""))

signaling.bor.fit <- train(Resistant ~ .,  data = signaling.bor.data, method = "glmnet",
                           trControl = fit.control, tuneLength = 5)
saveRDS(signaling.bor.fit, file = paste("models/", pre, "_live_signaling_bor_fit_run3.rds", sep = ""))

Bcl2.bor.fit <- train(Resistant ~ .,  data = Bcl2.bor.data, method = "glmnet",
                      trControl = fit.control, tuneLength = 5)
saveRDS(Bcl2.bor.fit, file = paste("models/", pre, "_live_Bcl2_bor_fit_run3.rds", sep = ""))

most.dex.fit <- train(Resistant ~ .,  data = most.dex.data, method = "glmnet",
                      trControl = fit.control, tuneLength = 5)
saveRDS(most.dex.fit, file = paste("models/", pre, "_live_most_dex_fit_run3.rds", sep = ""))

signaling.dex.fit <- train(Resistant ~ .,  data = signaling.dex.data, method = "glmnet",
                           trControl = fit.control, tuneLength = 5)
saveRDS(signaling.dex.fit, file = paste("models/", pre, "_live_signaling_dex_fit_run3.rds", sep = ""))

Bcl2.dex.fit <- train(Resistant ~ .,  data = Bcl2.dex.data, method = "glmnet",
                      trControl = fit.control, tuneLength = 5)
saveRDS(Bcl2.dex.fit, file = paste("models/", pre, "_live_Bcl2_dex_fit_run3.rds", sep = ""))

###############################################################
### Assess Models Using Signaling/Bcl-2 Family on Test Data ###
###############################################################

test.bor.data <- readRDS(paste("data/", pre, "_test_live_bortezomib_data_run3.rds", sep = ""))
test.dex.data <- readRDS(paste("data/", pre, "_test_live_dexamethasone_data_run3.rds", sep = ""))

# ROC Curve
most.bor.roc <- SaveROC(most.bor.fit, test.bor.data, response = "Resistant",
                        output.name = paste("roc_curves/", pre, "_live_most_bor_roc.rds", sep = ""),
                        ci = TRUE, smooth = FALSE)
most.dex.roc <- SaveROC(most.dex.fit, test.dex.data, response = "Resistant",
                        output.name = paste("roc_curves/", pre, "_live_most_dex_roc.rds", sep = ""),
                        ci = TRUE, smooth = FALSE)
signaling.bor.roc <- SaveROC(signaling.bor.fit, test.bor.data, response = "Resistant",
                             output.name = paste("roc_curves/", pre, "_live_signaling_bor_roc.rds", sep = ""),
                             ci = TRUE, smooth = FALSE)
signaling.dex.roc <- SaveROC(signaling.dex.fit, test.dex.data, response = "Resistant",
                             output.name = paste("roc_curves/", pre, "_live_signaling_dex_roc.rds", sep = ""),
                             ci = TRUE, smooth = FALSE)
Bcl2.bor.roc <- SaveROC(Bcl2.bor.fit, test.bor.data, response = "Resistant",
                        output.name = paste("roc_curves/", pre, "_live_Bcl2_bor_roc.rds", sep = ""),
                        ci = TRUE, smooth = FALSE)
Bcl2.dex.roc <- SaveROC(Bcl2.dex.fit, test.dex.data, response = "Resistant",
                        output.name = paste("roc_curves/", pre, "_live_Bcl2_dex_roc.rds", sep = ""),
                        ci = TRUE, smooth = FALSE)

###################################################
### Assess Dex/Bor Models on Opposite Drug Data ###
###################################################

# calculate test accuracy
# cross-test B model on D, and D model on B
# assess accuracy of classification

# ROC Curve
bor.on.dex.most.roc <- SaveROC(most.bor.fit, test.dex.data, response = "Resistant",
                               output.name = paste("roc_curves/", pre, "_live_most_bor_model_on_dex_data_roc.rds", sep = ""),
                               ci = TRUE, smooth = FALSE)
dex.on.bor.most.roc <- SaveROC(most.dex.fit, test.bor.data, response = "Resistant",
                               output.name = paste("roc_curves/", pre, "_live_most_dex_model_on_bor_data_roc.rds", sep = ""),
                               ci = TRUE, smooth = FALSE)

bor.on.dex.signaling.roc <- SaveROC(signaling.bor.fit, test.dex.data, response = "Resistant",
                                    output.name = paste("roc_curves/", pre, "_live_signaling_bor_model_on_dex_data_roc.rds", sep = ""),
                                    ci = TRUE, smooth = FALSE)
dex.on.bor.signaling.roc <- SaveROC(signaling.dex.fit, test.bor.data, response = "Resistant",
                                    output.name = paste("roc_curves/", pre, "_live_signaling_dex_model_on_bor_data_roc.rds", sep = ""),
                                    ci = TRUE, smooth = FALSE)

bor.on.dex.Bcl2.roc <- SaveROC(Bcl2.bor.fit, test.dex.data, response = "Resistant",
                               output.name = paste("roc_curves/", pre, "_live_Bcl2_bor_model_on_dex_data_roc.rds", sep = ""),
                               ci = TRUE, smooth = FALSE)
dex.on.bor.Bcl2.roc <- SaveROC(Bcl2.dex.fit, test.bor.data, response = "Resistant",
                               output.name = paste("roc_curves/", pre, "_live_Bcl2_dex_model_on_bor_data_roc.rds", sep = ""),
                               ci = TRUE, smooth = FALSE)

###############################################
### Determine aCasp3+/- Gates for DMSO Data ###
###############################################

all.files.together <- readRDS("/Users/mesako/Desktop/20180211_Bcl-2_Analysis/all_runs_live_cleaned.rds")
files.to.use <- names(all.files.together)[grepl(names(all.files.together), pattern = "WT")]
files.to.use <- files.to.use[grepl(files.to.use, pattern = "RunThree")]

DMSO.files <- c(files.to.use[grepl(files.to.use, pattern = "untreated")],
                files.to.use[grepl(files.to.use, pattern = "DMSO")])
DMSO.data <- all.files.together[DMSO.files]
rm(all.files.together)

DMSO.data <- CollapseBy(DMSO.data, c("00", "24", "72"),
                        "Timepoint", subsample = FALSE)

# Non-Normalized Values
bor.low.cutoff <- 2.25
bor.high.cutoff <- 3.25
dex.low.cutoff <- 2
dex.high.cutoff <- 3
xlim <- c(-1, 7)

qplot(DMSO.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Bortezomib",
      xlim = xlim) +
  geom_vline(xintercept = bor.low.cutoff, col = "red") + 
  geom_vline(xintercept = bor.high.cutoff, col = "red") +
  geom_vline(xintercept = dex.low.cutoff, col = "blue") + 
  geom_vline(xintercept = dex.high.cutoff, col = "blue") 
ggsave("run3_DMSO_alltime_cCaspase3.png", width = 6, height = 4, dpi = 400)

temp.DMSO.data <- DMSO.data[DMSO.data[, "Timepoint"] == "72", ]

qplot(temp.DMSO.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Bortezomib",
      xlim = xlim) +
  geom_vline(xintercept = bor.low.cutoff, col = "red") + 
  geom_vline(xintercept = bor.high.cutoff, col = "red") +
  geom_vline(xintercept = dex.low.cutoff, col = "blue") + 
  geom_vline(xintercept = dex.high.cutoff, col = "blue") 
ggsave("run3_DMSO_lasttime_cCaspase3.png", width = 6, height = 4, dpi = 400)

################################################
### Assess Most/ScoreDiff Model on DMSO Data ###
################################################

all.files.together <- readRDS("/Users/mesako/Desktop/20180211_Bcl-2_Analysis/all_runs_live_cleaned.rds")
files.to.use <- names(all.files.together)[grepl(names(all.files.together), pattern = "WT")]
files.to.use <- files.to.use[grepl(files.to.use, pattern = "RunThree")]

if (analysis.mode == "alltime") {
  DMSO.files <- c(files.to.use[grepl(files.to.use, pattern = "untreated")],
                  files.to.use[grepl(files.to.use, pattern = "DMSO")])
  DMSO.data <- all.files.together[DMSO.files]
  DMSO.data <- CollapseBy(DMSO.data, c("00", "24", "72"),
                          "Timepoint", subsample = FALSE)
} else if (analysis.mode == "lasttime") {
  files.to.use <- files.to.use[grepl(files.to.use, pattern = "DMSO-")]
  DMSO.files <- files.to.use[grepl(files.to.use, pattern = "72")]
  DMSO.data <- all.files.together[[DMSO.files]]
} else {
  stop("mode not recognized!")
}
rm(all.files.together)

all.files.together.norm <- readRDS("/Users/mesako/Desktop/20180211_Bcl-2_Analysis/all_runs_live_cleaned_quant_norm.rds")
if (analysis.mode == "alltime") {
  DMSO.files <- c(files.to.use[grepl(files.to.use, pattern = "untreated")],
                  files.to.use[grepl(files.to.use, pattern = "DMSO")])
  norm.DMSO.data <- all.files.together.norm[DMSO.files]
  norm.DMSO.data <- CollapseBy(norm.DMSO.data, c("00", "24", "72"),
                               "Timepoint", subsample = FALSE)
} else if (analysis.mode == "lasttime") {
  files.to.use <- files.to.use[grepl(files.to.use, pattern = "DMSO-")]
  DMSO.files <- files.to.use[grepl(files.to.use, pattern = "72")]
  norm.DMSO.data <- all.files.together.norm[[DMSO.files]]
} else {
  stop("mode not recognized!")
}
rm(all.files.together.norm)

# add Resistant TRUE/FALSE status to all cells
which.ind <- which(DMSO.data$aCasp3 < min(bor.low.cutoff, dex.low.cutoff))
nonapoptotic.DMSO.data <- norm.DMSO.data[which.ind, ]
which.ind <- which(DMSO.data$aCasp3 > max(bor.high.cutoff, dex.high.cutoff))
apoptotic.DMSO.data <- norm.DMSO.data[which.ind, ]
nonapoptotic.DMSO.data <- cbind(nonapoptotic.DMSO.data, Status = rep("Resistant", times = nrow(nonapoptotic.DMSO.data)))
apoptotic.DMSO.data <- cbind(apoptotic.DMSO.data, Status = rep("Sensitive", times = nrow(apoptotic.DMSO.data)))
recombined.DMSO.data <- rbind(nonapoptotic.DMSO.data, apoptotic.DMSO.data)
recombined.DMSO.data <- recombined.DMSO.data[sample(1:nrow(recombined.DMSO.data), size = nrow(recombined.DMSO.data)), ]

recombined.DMSO.data <- ConvertClassto01(recombined.DMSO.data, column = "Status", true.label = "Resistant")
recombined.DMSO.data[, "Resistant"] <- as.factor(recombined.DMSO.data[, "Resistant"])

setwd(paste(pre.folder, "run3_live_classification_models", sep = ""))
setwd(paste("run3_live_", pre,"_classification_models", sep = ""))
saveRDS(recombined.DMSO.data, paste("data/", pre, "_run3_DMSO_data.rds", sep = ""))

#################################################
### Assess Classification Models on DMSO Data ###
#################################################

most.bor.fit <- readRDS(paste("models/", pre, "_live_most_bor_fit_run3.rds", sep = ""))
most.dex.fit <- readRDS(paste("models/", pre, "_live_most_dex_fit_run3.rds", sep = ""))
score.diff.bor.fit <- readRDS(paste("models/", pre, "_live_scorediff_bor_fit_run3.rds", sep = ""))
score.diff.dex.fit <- readRDS(paste("models/", pre, "_live_scorediff_dex_fit_run3.rds", sep = ""))

# ROC Curve
most.bor.on.DMSO.roc <- SaveROC(most.bor.fit, recombined.DMSO.data, response = "Resistant",
                                output.name = paste("roc_curves/", pre, "_live_most_bor_on_DMSO_roc.rds", sep = ""),
                                ci = TRUE, smooth = FALSE)
most.dex.on.DMSO.roc <- SaveROC(most.dex.fit, recombined.DMSO.data, response = "Resistant",
                                output.name = paste("roc_curves/", pre, "_live_most_dex_on_DMSO_roc.rds", sep = ""),
                                ci = TRUE, smooth = FALSE)

most.bor.roc <- readRDS(paste("roc_curves/", pre, "_live_most_bor_roc.rds", sep = ""))
most.dex.roc <- readRDS(paste("roc_curves/", pre, "_live_most_dex_roc.rds", sep = ""))
most.bor.on.DMSO.roc <- readRDS(paste("roc_curves/", pre, "_live_most_bor_on_DMSO_roc.rds", sep = ""))
most.dex.on.DMSO.roc <- readRDS(paste("roc_curves/", pre, "_live_most_dex_on_DMSO_roc.rds", sep = ""))

########################################################################
### Assess Dex Most/ScoreDiff Model on Dex BimKO and BaxBak DKO Data ###
########################################################################

all.files.together.norm <- readRDS("/Users/mesako/Desktop/20180211_Bcl-2_Analysis/all_runs_live_cleaned_quant_norm.rds")
files.to.use <- names(all.files.together.norm)[grepl(names(all.files.together.norm), pattern = "RunThree")]
files.to.use <- files.to.use[grepl(files.to.use, pattern = "KO")]
files.to.use1 <- files.to.use[grepl(files.to.use, pattern = "_D-")]
files.to.use2 <- files.to.use[grepl(files.to.use, pattern = "_untreated-")]
files.to.use <- c(files.to.use1, files.to.use2)

if (analysis.mode == "alltime") {
  BakBaxDKO.files <- files.to.use[grepl(files.to.use, pattern = "Bak")]
  norm.BakBaxDKO.data <- all.files.together.norm[BakBaxDKO.files]
  norm.BakBaxDKO.data <- CollapseBy(norm.BakBaxDKO.data, c("00", "72"),
                                    "Timepoint", subsample = FALSE)
} else if (analysis.mode == "lasttime") {
  files.to.use <- files.to.use[grepl(files.to.use, pattern = "72")]
  norm.BakBaxDKO.data <- all.files.together.norm[[files.to.use[grepl(files.to.use, pattern = "Bak")]]]
} else {
  stop("mode not recognized!")
}
rm(all.files.together.norm)

all.files.together <- readRDS("/Users/mesako/Desktop/20180211_Bcl-2_Analysis/all_runs_live_cleaned.rds")
if (analysis.mode == "alltime") {
  BakBaxDKO.data <- all.files.together[BakBaxDKO.files]
  BakBaxDKO.data <- CollapseBy(BakBaxDKO.data, c("00", "72"),
                               "Timepoint", subsample = FALSE)
} else if (analysis.mode == "lasttime") {
  BakBaxDKO.data <- all.files.together[[files.to.use[grepl(files.to.use, pattern = "Bak")]]]
} else {
  stop("mode not recognized!")
}
rm(all.files.together)

# add Resistant TRUE/FALSE status to all cells
which.ind <- which(BakBaxDKO.data$aCasp3 < dex.low.cutoff)
nonapoptotic.BakBaxDKO.data <- norm.BakBaxDKO.data[which.ind, ]
which.ind <- which(BakBaxDKO.data$aCasp3 > dex.high.cutoff)
apoptotic.BakBaxDKO.data <- norm.BakBaxDKO.data[which.ind, ]
nonapoptotic.BakBaxDKO.data <- cbind(nonapoptotic.BakBaxDKO.data, Status = rep("Resistant", times = nrow(nonapoptotic.BakBaxDKO.data)))
apoptotic.BakBaxDKO.data <- cbind(apoptotic.BakBaxDKO.data, Status = rep("Sensitive", times = nrow(apoptotic.BakBaxDKO.data)))
recombined.BakBaxDKO.data <- rbind(nonapoptotic.BakBaxDKO.data, apoptotic.BakBaxDKO.data)
recombined.BakBaxDKO.data <- recombined.BakBaxDKO.data[sample(1:nrow(recombined.BakBaxDKO.data), size = nrow(recombined.BakBaxDKO.data)), ]

recombined.BakBaxDKO.data <- ConvertClassto01(recombined.BakBaxDKO.data, column = "Status", true.label = "Resistant")

recombined.BakBaxDKO.data[, "Resistant"] <- as.factor(recombined.BakBaxDKO.data[, "Resistant"])

saveRDS(recombined.BakBaxDKO.data, paste("data/", pre, "_run3_BakBaxDKO_data.rds", sep = ""))

################################################################
### Determine aCasp3+/- Gates for Other Run Data (Run 1 + 2) ###
################################################################

all.files.together <- readRDS("/Users/mesako/Desktop/20180211_Bcl-2_Analysis/all_runs_live_cleaned.rds")
files.to.use <- names(all.files.together)[grepl(names(all.files.together), pattern = "WT")]
files1.to.use <- files.to.use[grepl(files.to.use, pattern = "RunOne")]
files2.to.use <- files.to.use[grepl(files.to.use, pattern = "RunTwo")]

run1.bor.files <- c(files1.to.use[grepl(files1.to.use, pattern = "untreated")],
                    files1.to.use[grepl(files1.to.use, pattern = "B-")])
run2.bor.files <- c(files2.to.use[grepl(files2.to.use, pattern = "untreated")],
                    files2.to.use[grepl(files2.to.use, pattern = "B-")])
run1.dex.files <- c(files1.to.use[grepl(files1.to.use, pattern = "untreated")],
                    files1.to.use[grepl(files1.to.use, pattern = "D-")])
run2.dex.files <- c(files2.to.use[grepl(files2.to.use, pattern = "untreated")],
                    files2.to.use[grepl(files2.to.use, pattern = "D-")])

run1.bor.data <- all.files.together[run1.bor.files]
run2.bor.data <- all.files.together[run2.bor.files]
run1.dex.data <- all.files.together[run1.dex.files]
run2.dex.data <- all.files.together[run2.dex.files]
rm(all.files.together)

run1.bor.data <- CollapseBy(run1.bor.data, c("00", "06", "24"),
                            "Timepoint", subsample = FALSE)
run2.bor.data <- CollapseBy(run2.bor.data, c("00", "06", "24"),
                            "Timepoint", subsample = FALSE)
run1.dex.data <- CollapseBy(run1.dex.data, c("00", "24", "72"),
                            "Timepoint", subsample = FALSE)
run2.dex.data <- CollapseBy(run2.dex.data, c("00", "24", "72"),
                            "Timepoint", subsample = FALSE)

qplot(run1.bor.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Bortezomib",
      xlim = xlim) +
  geom_vline(xintercept = bor.low.cutoff, col = "red") + 
  geom_vline(xintercept = bor.high.cutoff, col = "red") 
ggsave("run1_bortezomib_alltime_cCaspase3.png", width = 6, height = 4, dpi = 400)

qplot(run2.bor.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Bortezomib",
      xlim = xlim) +
  geom_vline(xintercept = bor.low.cutoff, col = "red") + 
  geom_vline(xintercept = bor.high.cutoff, col = "red") 
ggsave("run2_bortezomib_alltime_cCaspase3.png", width = 6, height = 4, dpi = 400)

qplot(run1.dex.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Dexamethasone",
      xlim = xlim) +
  geom_vline(xintercept = dex.low.cutoff, col = "red") + 
  geom_vline(xintercept = dex.high.cutoff, col = "red") 
ggsave("run1_dexamethasone_alltime_cCaspase3.png", width = 6, height = 4, dpi = 400)

qplot(run2.dex.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Dexamethasone",
      xlim = xlim) +
  geom_vline(xintercept = dex.low.cutoff, col = "red") + 
  geom_vline(xintercept = dex.high.cutoff, col = "red") 
ggsave("run2_dexamethasone_alltime_cCaspase3.png", width = 6, height = 4, dpi = 400)

temp1.bor.data <- run1.bor.data[run1.bor.data[, "Timepoint"] == "24", ]
temp2.bor.data <- run2.bor.data[run2.bor.data[, "Timepoint"] == "24", ]
temp1.dex.data <- run1.dex.data[run1.dex.data[, "Timepoint"] == "72", ]
temp2.dex.data <- run2.dex.data[run2.dex.data[, "Timepoint"] == "72", ]

qplot(temp1.bor.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Bortezomib",
      xlim = xlim) +
  geom_vline(xintercept = bor.low.cutoff, col = "red") + 
  geom_vline(xintercept = bor.high.cutoff, col = "red") 
ggsave("run1_bortezomib_lasttime_cCaspase3.png", width = 6, height = 4, dpi = 400)

qplot(temp2.bor.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Bortezomib",
      xlim = xlim) +
  geom_vline(xintercept = bor.low.cutoff, col = "red") + 
  geom_vline(xintercept = bor.high.cutoff, col = "red") 
ggsave("run2_bortezomib_lasttime_cCaspase3.png", width = 6, height = 4, dpi = 400)

qplot(temp1.dex.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Dexamethasone",
      xlim = xlim) +
  geom_vline(xintercept = dex.low.cutoff, col = "red") + 
  geom_vline(xintercept = dex.high.cutoff, col = "red") 
ggsave("run1_dexamethasone_lasttime_cCaspase3.png", width = 6, height = 4, dpi = 400)

qplot(temp2.dex.data[, "aCasp3"], geom = "histogram",
      fill = I("grey"), col = I("black"),
      binwidth = 0.1, xlab = "aCasp3", main = "Dexamethasone",
      xlim = xlim) +
  geom_vline(xintercept = dex.low.cutoff, col = "red") + 
  geom_vline(xintercept = dex.high.cutoff, col = "red") 
ggsave("run2_dexamethasone_lasttime_cCaspase3.png", width = 6, height = 4, dpi = 400)

#################################################################
### Assess Most/ScoreDiff Model on Other Run Data (Run 1 + 2) ###
#################################################################

all.files.together <- readRDS("/Users/mesako/Desktop/20180211_Bcl-2_Analysis/all_runs_live_cleaned.rds")
files.to.use <- names(all.files.together)[grepl(names(all.files.together), pattern = "WT")]
files1.to.use <- files.to.use[grepl(files.to.use, pattern = "RunOne")]
files2.to.use <- files.to.use[grepl(files.to.use, pattern = "RunTwo")]

run1.bor.files <- c(files1.to.use[grepl(files1.to.use, pattern = "untreated")],
                    files1.to.use[grepl(files1.to.use, pattern = "B-")])
run2.bor.files <- c(files2.to.use[grepl(files2.to.use, pattern = "untreated")],
                    files2.to.use[grepl(files2.to.use, pattern = "B-")])
run1.dex.files <- c(files1.to.use[grepl(files1.to.use, pattern = "untreated")],
                    files1.to.use[grepl(files1.to.use, pattern = "D-")])
run2.dex.files <- c(files2.to.use[grepl(files2.to.use, pattern = "untreated")],
                    files2.to.use[grepl(files2.to.use, pattern = "D-")])

if (analysis.mode == "alltime") {
  run1.bor.data <- all.files.together[run1.bor.files]
  run2.bor.data <- all.files.together[run2.bor.files]
  run1.dex.data <- all.files.together[run1.dex.files]
  run2.dex.data <- all.files.together[run2.dex.files]
  run1.bor.data <- CollapseBy(run1.bor.data, c("00", "02", "03"),
                              "Timepoint", subsample = FALSE)
  run2.bor.data <- CollapseBy(run2.bor.data, c("00", "02", "03"),
                              "Timepoint", subsample = FALSE)
  run1.dex.data <- CollapseBy(run1.dex.data, c("00", "02", "03"),
                              "Timepoint", subsample = FALSE)
  run2.dex.data <- CollapseBy(run2.dex.data, c("00", "02", "03"),
                              "Timepoint", subsample = FALSE)
} else if (analysis.mode == "lasttime") {
  run1.bor.data <- all.files.together[[run1.bor.files[grepl(run1.bor.files, pattern = "24")]]]
  run2.bor.data <- all.files.together[[run2.bor.files[grepl(run2.bor.files, pattern = "24")]]]
  run1.dex.data <- all.files.together[[run1.dex.files[grepl(run1.dex.files, pattern = "72")]]]
  run2.dex.data <- all.files.together[[run2.dex.files[grepl(run2.dex.files, pattern = "72")]]]
} else {
  stop("mode not recognized!")
}
rm(all.files.together)

all.files.together.norm <- readRDS("/Users/mesako/Desktop/20180211_Bcl-2_Analysis/all_runs_live_cleaned_quant_norm.rds")
if (analysis.mode == "alltime") {
  norm.run1.bor.data <- all.files.together.norm[run1.bor.files]
  norm.run2.bor.data <- all.files.together.norm[run2.bor.files]
  norm.run1.dex.data <- all.files.together.norm[run1.dex.files]
  norm.run2.dex.data <- all.files.together.norm[run2.dex.files]
  norm.run1.bor.data <- CollapseBy(norm.run1.bor.data, c("00", "02", "03"),
                                   "Timepoint", subsample = FALSE)
  norm.run2.bor.data <- CollapseBy(norm.run2.bor.data, c("00", "02", "03"),
                                   "Timepoint", subsample = FALSE)
  norm.run1.dex.data <- CollapseBy(norm.run1.dex.data, c("00", "02", "03"),
                                   "Timepoint", subsample = FALSE)
  norm.run2.dex.data <- CollapseBy(norm.run2.dex.data, c("00", "02", "03"),
                                   "Timepoint", subsample = FALSE)
} else if (analysis.mode == "lasttime") {
  temp.run1.bor.files <- run1.bor.files[grepl(run1.bor.files, pattern = "24")]
  temp.run2.bor.files <- run2.bor.files[grepl(run2.bor.files, pattern = "24")]
  temp.run1.dex.files <- run1.dex.files[grepl(run1.dex.files, pattern = "72")]
  temp.run2.dex.files <- run2.dex.files[grepl(run2.dex.files, pattern = "72")]
  norm.run1.bor.data <- all.files.together.norm[[temp.run1.bor.files]]
  norm.run2.bor.data <- all.files.together.norm[[temp.run2.bor.files]]
  norm.run1.dex.data <- all.files.together.norm[[temp.run1.dex.files]]
  norm.run2.dex.data <- all.files.together.norm[[temp.run2.dex.files]]
} else {
  stop("mode not recognized!")
}
rm(all.files.together.norm)

# add Resistant TRUE/FALSE status to all cells
which.ind <- which(run1.bor.data$aCasp3 < bor.low.cutoff)
nonapoptotic.run1.bor.data <- norm.run1.bor.data[which.ind, ]
which.ind <- which(run1.bor.data$aCasp3 > bor.high.cutoff)
apoptotic.run1.bor.data <- norm.run1.bor.data[which.ind, ]
nonapoptotic.run1.bor.data <- cbind(nonapoptotic.run1.bor.data, Status = rep("Resistant", times = nrow(nonapoptotic.run1.bor.data)))
apoptotic.run1.bor.data <- cbind(apoptotic.run1.bor.data, Status = rep("Sensitive", times = nrow(apoptotic.run1.bor.data)))
recombined.run1.bor.data <- rbind(nonapoptotic.run1.bor.data, apoptotic.run1.bor.data)
recombined.run1.bor.data <- recombined.run1.bor.data[sample(1:nrow(recombined.run1.bor.data), size = nrow(recombined.run1.bor.data)), ]

which.ind <- which(run2.bor.data$aCasp3 < bor.low.cutoff)
nonapoptotic.run2.bor.data <- norm.run2.bor.data[which.ind, ]
which.ind <- which(run2.bor.data$aCasp3 > bor.high.cutoff)
apoptotic.run2.bor.data <- norm.run2.bor.data[which.ind, ]
nonapoptotic.run2.bor.data <- cbind(nonapoptotic.run2.bor.data, Status = rep("Resistant", times = nrow(nonapoptotic.run2.bor.data)))
apoptotic.run2.bor.data <- cbind(apoptotic.run2.bor.data, Status = rep("Sensitive", times = nrow(apoptotic.run2.bor.data)))
recombined.run2.bor.data <- rbind(nonapoptotic.run2.bor.data, apoptotic.run2.bor.data)
recombined.run2.bor.data <- recombined.run2.bor.data[sample(1:nrow(recombined.run2.bor.data), size = nrow(recombined.run2.bor.data)), ]

which.ind <- which(run1.dex.data$aCasp3 < dex.low.cutoff)
nonapoptotic.run1.dex.data <- norm.run1.dex.data[which.ind, ]
which.ind <- which(run1.dex.data$aCasp3 > dex.high.cutoff)
apoptotic.run1.dex.data <- norm.run1.dex.data[which.ind, ]
nonapoptotic.run1.dex.data <- cbind(nonapoptotic.run1.dex.data, Status = rep("Resistant", times = nrow(nonapoptotic.run1.dex.data)))
apoptotic.run1.dex.data <- cbind(apoptotic.run1.dex.data, Status = rep("Sensitive", times = nrow(apoptotic.run1.dex.data)))
recombined.run1.dex.data <- rbind(nonapoptotic.run1.dex.data, apoptotic.run1.dex.data)
recombined.run1.dex.data <- recombined.run1.dex.data[sample(1:nrow(recombined.run1.dex.data), size = nrow(recombined.run1.dex.data)), ]

which.ind <- which(run2.dex.data$aCasp3 < dex.low.cutoff)
nonapoptotic.run2.dex.data <- norm.run2.dex.data[which.ind, ]
which.ind <- which(run2.dex.data$aCasp3 > dex.high.cutoff)
apoptotic.run2.dex.data <- norm.run2.dex.data[which.ind, ]
nonapoptotic.run2.dex.data <- cbind(nonapoptotic.run2.dex.data, Status = rep("Resistant", times = nrow(nonapoptotic.run2.dex.data)))
apoptotic.run2.dex.data <- cbind(apoptotic.run2.dex.data, Status = rep("Sensitive", times = nrow(apoptotic.run2.dex.data)))
recombined.run2.dex.data <- rbind(nonapoptotic.run2.dex.data, apoptotic.run2.dex.data)
recombined.run2.dex.data <- recombined.run2.dex.data[sample(1:nrow(recombined.run2.dex.data), size = nrow(recombined.run2.dex.data)), ]

recombined.run1.bor.data <- ConvertClassto01(recombined.run1.bor.data, column = "Status", true.label = "Resistant")
recombined.run2.bor.data <- ConvertClassto01(recombined.run2.bor.data, column = "Status", true.label = "Resistant")
recombined.run1.dex.data <- ConvertClassto01(recombined.run1.dex.data, column = "Status", true.label = "Resistant")
recombined.run2.dex.data <- ConvertClassto01(recombined.run2.dex.data, column = "Status", true.label = "Resistant")

recombined.run1.bor.data[, "Resistant"] <- as.factor(recombined.run1.bor.data[, "Resistant"])
recombined.run2.bor.data[, "Resistant"] <- as.factor(recombined.run2.bor.data[, "Resistant"])
recombined.run1.dex.data[, "Resistant"] <- as.factor(recombined.run1.dex.data[, "Resistant"])
recombined.run2.dex.data[, "Resistant"] <- as.factor(recombined.run2.dex.data[, "Resistant"])

saveRDS(recombined.run1.bor.data, paste("data/", pre, "_run1_bortezomib_data.rds", sep = ""))
saveRDS(recombined.run2.bor.data, paste("data/", pre, "_run2_bortezomib_data.rds", sep = ""))
saveRDS(recombined.run1.dex.data, paste("data/", pre, "_run1_dexamethasone_data.rds", sep = ""))
saveRDS(recombined.run2.dex.data, paste("data/", pre, "_run2_dexamethasone_data.rds", sep = ""))

# recombined.run1.bor.data <- readRDS(paste("data/", pre, "_run1_bortezomib_data.rds", sep = ""))
# recombined.run2.bor.data <- readRDS(paste("data/", pre, "_run2_bortezomib_data.rds", sep = ""))
# recombined.run1.dex.data <- readRDS(paste("data/", pre, "_run1_dexamethasone_data.rds", sep = ""))
# recombined.run2.dex.data <- readRDS(paste("data/", pre, "_run2_dexamethasone_data.rds", sep = ""))

# most.bor.fit <- readRDS(paste("models/", pre, "_live_most_bor_fit_run3.rds", sep = ""))
# most.dex.fit <- readRDS(paste("models/", pre, "_live_most_dex_fit_run3.rds", sep = ""))
# score.diff.bor.fit <- readRDS(paste("models/", pre, "_live_scorediff_bor_fit_run3.rds", sep = ""))
# score.diff.dex.fit <- readRDS(paste("models/", pre, "_live_scorediff_dex_fit_run3.rds", sep = ""))

# ROC Curve
run1.most.bor.roc <- SaveROC(most.bor.fit, recombined.run1.bor.data, response = "Resistant",
                             output.name = paste("roc_curves/", pre, "_live_run1_most_bor_roc.rds", sep = ""),
                             ci = TRUE, smooth = FALSE)
run2.most.bor.roc <- SaveROC(most.bor.fit, recombined.run2.bor.data, response = "Resistant",
                             output.name = paste("roc_curves/", pre, "_live_run2_most_bor_roc.rds", sep = ""),
                             ci = TRUE, smooth = FALSE)
run1.most.dex.roc <- SaveROC(most.dex.fit, recombined.run1.dex.data, response = "Resistant",
                             output.name = paste("roc_curves/", pre, "_live_run1_most_dex_roc.rds", sep = ""),
                             ci = TRUE, smooth = FALSE)
run2.most.dex.roc <- SaveROC(most.dex.fit, recombined.run2.dex.data, response = "Resistant",
                             output.name = paste("roc_curves/", pre, "_live_run2_most_dex_roc.rds", sep = ""),
                             ci = TRUE, smooth = FALSE)

run3.most.bor.roc <- readRDS(paste("roc_curves/", pre, "_live_most_bor_roc.rds", sep = ""))
run3.most.dex.roc <- readRDS(paste("roc_curves/", pre, "_live_most_dex_roc.rds", sep = ""))

#############################
### Load Data (if needed) ###
#############################

setwd(paste(pre.folder, "run3_live_classification_models", sep = ""))
setwd(paste("run3_live_", pre,"_classification_models", sep = ""))
training.bor.data <- readRDS(paste("data/", pre, "_training_live_bortezomib_data_run3.rds", sep = ""))
training.dex.data <- readRDS(paste("data/", pre, "_training_live_dexamethasone_data_run3.rds", sep = ""))

signaling.channels <- c("Akt_pS473", "p53", "STAT5_pY694", "cMyc", "IkBa",
                        "S6_pS235", "CREB_pS133", "p38_pT180", "Erk_p44")
pro.survival.Bcl2 <- c("Mcl1", "Bcl2", "BclxL")
pro.apoptotic.Bcl2 <- c("Bim", "Bax", "Bak_activated")
Bcl2.family.channels <- c(pro.survival.Bcl2, pro.apoptotic.Bcl2)
most.channels <- c(Bcl2.family.channels, signaling.channels)

most.bor.data <- subset(training.bor.data, select = c(most.channels, "Resistant"))
most.dex.data <- subset(training.dex.data, select = c(most.channels, "Resistant"))
Bcl2.bor.data <- subset(training.bor.data, select = c(Bcl2.family.channels, "Resistant"))
Bcl2.dex.data <- subset(training.dex.data, select = c(Bcl2.family.channels, "Resistant"))

#########################################################
### Classification Tree Modeling For Resistant Status ###
#########################################################

# Random Forests
library(randomForest)

temp.fit.control <- trainControl(method = "repeatedcv", number = 10, repeats = 10)

Bcl2.rf.bor.fit <- train(Resistant ~ .,  data = Bcl2.bor.data, method = "rf",
                         trControl = temp.fit.control, ntree = 50)
saveRDS(Bcl2.rf.bor.fit, file = paste("models/", pre, "_live_Bcl2_bor_rf_fit_run3.rds", sep = ""))

Bcl2.rf.dex.fit <- train(Resistant ~ .,  data = Bcl2.dex.data, method = "rf",
                         trControl = temp.fit.control, ntree = 50)
saveRDS(Bcl2.rf.dex.fit, file = paste("models/", pre, "_live_Bcl2_dex_rf_fit_run3.rds", sep = ""))

most.rf.bor.fit <- train(Resistant ~ .,  data = most.bor.data, method = "rf",
                         trControl = temp.fit.control, ntree = 50)
saveRDS(most.rf.bor.fit, file = paste("models/", pre, "_live_most_bor_rf_fit_run3.rds", sep = ""))

most.rf.dex.fit <- train(Resistant ~ .,  data = most.dex.data, method = "rf",
                         trControl = temp.fit.control, ntree = 50)
saveRDS(most.rf.dex.fit, file = paste("models/", pre, "_live_most_dex_rf_fit_run3.rds", sep = ""))

######################################################
### Assess Classification Tree Models on Test Data ###
######################################################

test.bor.data <- readRDS(paste("data/", pre, "_test_live_bortezomib_data_run3.rds", sep = ""))
test.dex.data <- readRDS(paste("data/", pre, "_test_live_dexamethasone_data_run3.rds", sep = ""))

# ROC Curve
Bcl2.rf.bor.roc <- SaveROC(Bcl2.rf.bor.fit, test.bor.data, response = "Resistant",
                           output.name = paste("roc_curves/", pre, "_live_Bcl2_rf_bor_roc.rds", sep = ""),
                           ci = TRUE, smooth = FALSE)
Bcl2.rf.dex.roc <- SaveROC(Bcl2.rf.dex.fit, test.dex.data, response = "Resistant",
                           output.name = paste("roc_curves/", pre, "_live_Bcl2_rf_dex_roc.rds", sep = ""),
                           ci = TRUE, smooth = FALSE)
most.rf.bor.roc <- SaveROC(most.rf.bor.fit, test.bor.data, response = "Resistant",
                           output.name = paste("roc_curves/", pre, "_live_most_rf_bor_roc.rds", sep = ""),
                           ci = TRUE, smooth = FALSE)
most.rf.dex.roc <- SaveROC(most.rf.dex.fit, test.dex.data, response = "Resistant",
                           output.name = paste("roc_curves/", pre, "_live_most_rf_dex_roc.rds", sep = ""),
                           ci = TRUE, smooth = FALSE)

###############################################################
### Assess Classification Tree Models on Opposite Drug Data ###
###############################################################

test.bor.data <- readRDS(paste("data/", pre, "_test_live_bortezomib_data_run3.rds", sep = ""))
test.dex.data <- readRDS(paste("data/", pre, "_test_live_dexamethasone_data_run3.rds", sep = ""))

# ROC Curve
Bcl2.rf.bor.on.dex.roc <- SaveROC(Bcl2.rf.bor.fit, test.dex.data, response = "Resistant",
                                  output.name = paste("roc_curves/", pre, "_live_Bcl2_rf_bor_on_dex_roc.rds", sep = ""),
                                  ci = TRUE, smooth = FALSE)
Bcl2.rf.dex.on.bor.roc <- SaveROC(Bcl2.rf.dex.fit, test.bor.data, response = "Resistant",
                                  output.name = paste("roc_curves/", pre, "_live_Bcl2_rf_dex_on_bor_roc.rds", sep = ""),
                                  ci = TRUE, smooth = FALSE)
most.rf.bor.on.dex.roc <- SaveROC(most.rf.bor.fit, test.dex.data, response = "Resistant",
                                  output.name = paste("roc_curves/", pre, "_live_most_rf_bor_on_dex_roc.rds", sep = ""),
                                  ci = TRUE, smooth = FALSE)
most.rf.dex.on.bor.roc <- SaveROC(most.rf.dex.fit, test.bor.data, response = "Resistant",
                                  output.name = paste("roc_curves/", pre, "_live_most_rf_dex_on_bor_roc.rds", sep = ""),
                                  ci = TRUE, smooth = FALSE)

#######################################################################
### Assess Classification Tree Models on Other Run Data (Run 1 + 2) ###
#######################################################################

recombined.run1.bor.data <- readRDS(paste("data/", pre, "_run1_bortezomib_data.rds", sep = ""))
recombined.run2.bor.data <- readRDS(paste("data/", pre, "_run2_bortezomib_data.rds", sep = ""))
recombined.run1.dex.data <- readRDS(paste("data/", pre, "_run1_dexamethasone_data.rds", sep = ""))
recombined.run2.dex.data <- readRDS(paste("data/", pre, "_run2_dexamethasone_data.rds", sep = ""))

# ROC Curve
run1.rf.Bcl2.bor.roc <- SaveROC(Bcl2.rf.bor.fit, recombined.run1.bor.data, response = "Resistant",
                                output.name = paste("roc_curves/", pre, "_live_run1_rf_Bcl2_bor_roc.rds", sep = ""),
                                ci = TRUE, smooth = FALSE)
run2.rf.Bcl2.bor.roc <- SaveROC(Bcl2.rf.bor.fit, recombined.run2.bor.data, response = "Resistant",
                                output.name = paste("roc_curves/", pre, "_live_run2_rf_Bcl2_bor_roc.rds", sep = ""),
                                ci = TRUE, smooth = FALSE)
run1.rf.Bcl2.dex.roc <- SaveROC(Bcl2.rf.dex.fit, recombined.run1.dex.data, response = "Resistant",
                                output.name = paste("roc_curves/", pre, "_live_run1_rf_Bcl2_dex_roc.rds", sep = ""),
                                ci = TRUE, smooth = FALSE)
run2.rf.Bcl2.dex.roc <- SaveROC(Bcl2.rf.dex.fit, recombined.run2.dex.data, response = "Resistant",
                                output.name = paste("roc_curves/", pre, "_live_run2_rf_Bcl2_dex_roc.rds", sep = ""),
                                ci = TRUE, smooth = FALSE)

run1.rf.most.bor.roc <- SaveROC(most.rf.bor.fit, recombined.run1.bor.data, response = "Resistant",
                                output.name = paste("roc_curves/", pre, "_live_run1_rf_most_bor_roc.rds", sep = ""),
                                ci = TRUE, smooth = FALSE)
run2.rf.most.bor.roc <- SaveROC(most.rf.bor.fit, recombined.run2.bor.data, response = "Resistant",
                                output.name = paste("roc_curves/", pre, "_live_run2_rf_most_bor_roc.rds", sep = ""),
                                ci = TRUE, smooth = FALSE)
run1.rf.most.dex.roc <- SaveROC(most.rf.dex.fit, recombined.run1.dex.data, response = "Resistant",
                                output.name = paste("roc_curves/", pre, "_live_run1_rf_most_dex_roc.rds", sep = ""),
                                ci = TRUE, smooth = FALSE)
run2.rf.most.dex.roc <- SaveROC(most.rf.dex.fit, recombined.run2.dex.data, response = "Resistant",
                                output.name = paste("roc_curves/", pre, "_live_run2_rf_most_dex_roc.rds", sep = ""),
                                ci = TRUE, smooth = FALSE)

######################################################
### Assess Classification Tree Models on DMSO Data ###
######################################################

recombined.DMSO.data <- readRDS(paste("data/", pre, "_run3_DMSO_data.rds", sep = ""))

Bcl2.rf.bor.fit <- readRDS(paste("models/", pre, "_live_Bcl2_bor_rf_fit_run3.rds", sep = ""))
Bcl2.rf.dex.fit <- readRDS(paste("models/", pre, "_live_Bcl2_dex_rf_fit_run3.rds", sep = ""))
most.rf.bor.fit <- readRDS(paste("models/", pre, "_live_most_bor_rf_fit_run3.rds", sep = ""))
most.rf.dex.fit <- readRDS(paste("models/", pre, "_live_most_dex_rf_fit_run3.rds", sep = ""))

# ROC Curve
Bcl2.rf.bor.on.DMSO.roc <- SaveROC(Bcl2.rf.bor.fit, recombined.DMSO.data, response = "Resistant",
                                   output.name = paste("roc_curves/", pre, "_live_Bcl2_rf_bor_on_DMSO_roc.rds", sep = ""),
                                   ci = TRUE, smooth = FALSE)
Bcl2.rf.dex.on.DMSO.roc <- SaveROC(Bcl2.rf.dex.fit, recombined.DMSO.data, response = "Resistant",
                                   output.name = paste("roc_curves/", pre, "_live_Bcl2_rf_dex_on_DMSO_roc.rds", sep = ""),
                                   ci = TRUE, smooth = FALSE)

most.rf.bor.on.DMSO.roc <- SaveROC(most.rf.bor.fit, recombined.DMSO.data, response = "Resistant",
                                   output.name = paste("roc_curves/", pre, "_live_most_rf_bor_on_DMSO_roc.rds", sep = ""),
                                   ci = TRUE, smooth = FALSE)
most.rf.dex.on.DMSO.roc <- SaveROC(most.rf.dex.fit, recombined.DMSO.data, response = "Resistant",
                                   output.name = paste("roc_curves/", pre, "_live_most_rf_dex_on_DMSO_roc.rds", sep = ""),
                                   ci = TRUE, smooth = FALSE)

#############################################
### Make Additional Random Forest Figures ###
#############################################

setwd(paste(pre.folder, "run3_live_classification_models", sep = ""))
setwd(paste("run3_live_", pre,"_classification_models", sep = ""))

Bcl2.rf.bor.fit <- readRDS(file = paste("models/", pre, "_live_Bcl2_bor_rf_fit_run3.rds", sep = ""))
Bcl2.rf.dex.fit <- readRDS(file = paste("models/", pre, "_live_Bcl2_dex_rf_fit_run3.rds", sep = ""))
most.rf.bor.fit <- readRDS(file = paste("models/", pre, "_live_most_bor_rf_fit_run3.rds", sep = ""))
most.rf.dex.fit <- readRDS(file = paste("models/", pre, "_live_most_dex_rf_fit_run3.rds", sep = ""))

file.name <-  paste(pre, "_Bcl2_rf_bor_fit_oob_vs_numtrees_run3.png", sep = "")
png(file.name, width = 960, height = 720, units = "px", res = 200)
plot(Bcl2.rf.bor.fit$finalModel, main = "Bcl2_rf_bor_fit")
dev.off()

file.name <-  paste(pre, "_Bcl2_rf_dex_fit_oob_vs_numtrees_run3.png", sep = "")
png(file.name, width = 960, height = 720, units = "px", res = 200)
plot(Bcl2.rf.dex.fit$finalModel, main = "Bcl2_rf_dex_fit")
dev.off()

file.name <-  paste(pre, "_most_rf_bor_fit_oob_vs_numtrees_run3.png", sep = "")
png(file.name, width = 960, height = 720, units = "px", res = 200)
plot(most.rf.bor.fit$finalModel, main = "most_rf_bor_fit")
dev.off()

file.name <-  paste(pre, "_most_rf_dex_fit_oob_vs_numtrees_run3.png", sep = "")
png(file.name, width = 960, height = 720, units = "px", res = 200)
plot(most.rf.dex.fit$finalModel, main = "most_rf_dex_fit")
dev.off()

most.rf.dex.fit <- readRDS(paste("models/", pre, "_live_most_dex_rf_fit_run3.rds", sep = ""))

recombined.BakBaxDKO.data <- readRDS(paste("data/", pre, "_run3_BakBaxDKO_data.rds", sep = ""))

most.rf.dex.on.BakBaxDKO.roc <- SaveROC(most.rf.dex.fit, recombined.BakBaxDKO.data, response = "Resistant",
                                        output.name = paste("roc_curves/", pre, "_live_accuracy_most_rf_dex_fit_on_BakBaxDKO_roc.rds", sep = ""),
                                        ci = TRUE, smooth = FALSE)

Bcl2.rf.dex.fit <- readRDS(paste("models/", pre, "_live_Bcl2_dex_rf_fit_run3.rds", sep = ""))

Bcl2.rf.dex.on.BakBaxDKO.roc <- SaveROC(Bcl2.rf.dex.fit, recombined.BakBaxDKO.data, response = "Resistant",
                                        output.name = paste("roc_curves/", pre, "_live_accuracy_Bcl2_rf_dex_fit_on_BakBaxDKO_roc.rds", sep = ""),
                                        ci = TRUE, smooth = FALSE)

rm(list = ls())

library(igraph)
library(flowCore)
library(robustbase)
library(ggplot2)
library(preprocessCore)
library(FLOWMAPR)

pre.folder <- "/Users/mesako/Desktop/20180211_Bcl-2_Analysis/"
setwd(pre.folder)

###################################################
### Load All FCS Files and Establish Parameters ###
###################################################

all.files.folder <- "/Users/mesako/Desktop/20180211_Bcl-2_Analysis/FCS_Files_Live_Together"
all.file.names <- FLOWMAPR::GetFCSNames(all.files.folder)
var.annotate <- FLOWMAPR::ConstructVarAnnotate(all.file.names[1])
var.remove <- FLOWMAPR::SuggestVarRemove(var.annotate)
var.remove <- c(var.remove, "File_Number", "Yb171Di", "CD61",
                "ZAP70_pY319", "Bi209Di")

############################################
### Organize All FCS Files For FLOW-MAPs ###
############################################

prefolder <- "/Users/mesako/Desktop/20180211_Bcl-2_Analysis/FCS_Files_for_Final_FLOWMAPs/"
all.files <- list()
all.modes <- c("single", "single", "single", "multi", "multi")

for (i in 1:3) {
  run.files.list <- list()
  # Singlets WT + Bortezomib singleFLOWMAP time course
  run.files.list[[1]] <- paste(prefolder, "Run", i, "/Singlets_WT_B", sep = "")
  # Singlets WT + Dexamethasone singleFLOWMAP time course
  run.files.list[[2]] <- paste(prefolder, "Run", i, "/Singlets_WT_D", sep = "")
  # Singlets WT + DMSO singleFLOWMAP time course
  run.files.list[[3]] <- paste(prefolder, "Run", i, "/Singlets_WT_DMSO", sep = "")
  # Live WT vs. BaxBakDKO + Dexamethasone multiFLOWMAP combined time course
  run.files.list[[4]] <- paste(prefolder, "Run", i, "/Live_WT_v_BakBaxDKO_D", sep = "")
  # Live WT + Dexamethasone vs. Bortezomib multiFLOWMAP combined time course
  run.files.list[[5]] <- paste(prefolder, "Run", i, "/Live_WT_B_v_D", sep = "")
  all.files[[i]] <- run.files.list
}

############################################
### Produce FLOW-MAPs for All Three Runs ###
############################################

save.folder <- "/Users/mesako/Desktop/20180211_Bcl-2_Analysis/20180312_Final_Bcl-2_FLOW-MAPs"
minimum <- 10
maximum <- 30
distance.metric <- "manhattan"
max.nodes <- 6000

pro.survival.Bcl2 <- c("Mcl1", "Bcl2", "BclxL")
pro.apoptotic.Bcl2 <- c("Bim", "Bax", "Bak_activated")

name.sort <- TRUE
downsample <- FALSE
savePDFs <- TRUE
which.palette <- "bluered"

seeds.to.use <- c(1, 2, 3)
all.runs <- c(1, 2, 3)

cat("start time:", format(Sys.time(), "%a %b %d %X %Y"), "\n")
for (n in all.runs) {
  for (x in 1:length(all.files[[n]])) {
    for (i in seeds.to.use) {
      cat("files are", all.files[[n]][[x]], "\n")
      num.files <- length(list.files(all.files[[n]][[x]], recursive = TRUE,
                                     pattern = "\\.fcs"))
      cluster.numbers <- max.nodes / num.files
      cat("cluster.numbers is", cluster.numbers, "\n")
      subsamples <- cluster.numbers * 5
      cat("subsamples is", subsamples, "\n")
      seed.X <- i
      cat("seed is", seed.X, "\n")
      mode <- all.modes[x]
      cat("mode is", mode, "\n")
      files <- all.files[[n]][[x]]
      if (grepl(basename(files), pattern = "Live")) {
        clustering.var <- c("cPARP", "p38_pT180", "Erk_p44", "Rb_pS807", "Akt_pS473", "BclxL",
                            "Bax", "Bak_activated", "Cyclin_B", "Bcl2", "STAT5_pY694",
                            "Mcl1", "cMyc", "IkBa", "Bim", "Cyclin_A",
                            "H3_pS28", "aCasp3", "p53", "S6_pS235", "CREB_pS133")
      } else {
        print("using Cisplatin in clustering.var")
        clustering.var <- c("cPARP", "p38_pT180", "Erk_p44", "Rb_pS807", "Akt_pS473", "BclxL",
                            "Bax", "Bak_activated", "Cyclin_B", "Bcl2", "STAT5_pY694",
                            "Mcl1", "cMyc", "IkBa", "Bim", "Cyclin_A",
                            "H3_pS28", "aCasp3", "p53", "S6_pS235", "CREB_pS133")
        clustering.var <- c(clustering.var, "Cisplatin")
      }
    }
  }
}
cat("end time:", format(Sys.time(), "%a %b %d %X %Y"), "\n")

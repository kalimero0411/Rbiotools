#!/usr/bin/env Rscript
packages = c("RColorBrewer","rgl","R.utils")
invisible(
  suppressMessages(
    sapply(packages,FUN = function(x) {
       if(!x %in% rownames(installed.packages())){
         cat("Installing package: ",x,"\n",sep = "")
         BiocManager::install(x,update = FALSE,ask = FALSE)
       }
      cat("#####   Loading package: ",x,"   #####\n",sep = "")
      library(x,character.only = TRUE)
    })))

options(stringsAsFactors = FALSE)
args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
setwd(args[["wd"]])
cat("Working directory: ",getwd(),"\n")
load("PCA_data.RData")
factor = args[["factor"]]
cat("Creating 3D PCA for factor ",factor,"\n",sep = "")
PCA_3D(PC_factor = factor)

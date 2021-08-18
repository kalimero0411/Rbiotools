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
if("list" %in% names(args)){
  cat("Available factors: ")
  cat(paste(colnames(experimental_design), collapse = ","),"\n")
}else{
if(! "factors" %in% names(args)){
  factors = colnames(experimental_design)
}else{
  factors = unlist(strsplit(args[["factors"]],split = ","))
}
for(factor in factors){
  cat("Creating 3D PCA for factor ",factor,"\n",sep = "")
  PCA_3D(PC_factor = factor)
}
}
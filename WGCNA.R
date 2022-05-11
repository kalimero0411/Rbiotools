##### SETUP #####
packages = c("purrr", "bindr", "Biobase", "BiocGenerics", "checkmate", "digest","dplyr", "dynamicTreeCut", "flashClust", 
             "ggplot2", "gplots", "ggpubr", "IRanges", "RColorBrewer", "robust", "tidyr", "WGCNA","parallel","R.utils","MatrixGenerics","pheatmap")

if(interactive()){
  packages = c(packages,"rChoiceDialogs","BiocParallel")
}

invisible(
  suppressMessages(
    sapply(packages,FUN = function(x) {
    #   if(!x %in% rownames(installed.packages())){
    #     cat("Installing package: ",x,"\n",sep = "")
    #     BiocManager::install(x,update = FALSE,ask = FALSE)
    #   }
      cat("#####   Loading package: ",x,"   #####\n",sep = "")
      library(x,character.only = TRUE)
    })))

options(stringsAsFactors = FALSE)

if(!interactive()){
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  must_args = c("wd","name","threads")
  if(!all(must_args %in% names(args)) & !any(c("rds","tsv","csv") %in% names(args))){
    help = matrix(data = c("--rds    ","RDS file path",
                           "--factors    ","Factors for experimental design",
                           "--tsv","Data in tab-seperated format",
                           "--csv","Data in comma-seperated format",
                           "--wd","Working directory path",
                           "--name","Experiment name",
                           "--power","Scale Free Topology Model Fit",
                           "--TOM","Specify sft value and create Topology overlap matrix (TOM)",
                           "--blockwise","Perform blockwise calculation for TOM",
                           "--minmod","Minimum module size",
                           "--maxblock","Maximum block size (only used in blockwise)",
                           "--min_exp","Minimum expression for all samples of a gene",
                           "--threads","Number of CPU threads",
                           "--arg","Additional R arguments (multiple arguments in separate flags)")
                  ,ncol = 2,byrow = TRUE)
    prmatrix(help,quote = FALSE,rowlab = rep("",nrow(help)),collab = rep("",2))
    stop(paste0("Missing command line input --> rds/tsv/csv | ",paste(must_args[!must_args %in% names(args)],collapse = " | ")), call. = TRUE)
  }
  
  # saveRDS(data,file = "data.rds")
  cat("Loading RData file: ",args[["rdata"]],"\n",sep = "")
  # load(args[["rdata"]])
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  section = c()
  wd = args[["wd"]]
  Experiment_name = args[["name"]]
  blockwise = "blockwise" %in% names(args)
  if("power" %in% names(args)){
    section = c(section,"power")
  }
  if("TOM" %in% names(args)){
    softPower = as.numeric(args[["TOM"]])
    section = c(section,"TOM")
  }
  if("minmod" %in% names(args)){
    minModuleSize = as.numeric(args[["minmod"]])
  }else{
    minModuleSize = 50
  }
  if("maxblock" %in% names(args)){
    maxBlockSize = args[["maxblock"]]
  }else{
    maxBlockSize = 10000
  }
  if("min_exp" %in% names(args)){
    min_expression = as.numeric(args[["min_exp"]])
  }else{
    min_expression = 0
  }
  if("factors" %in% names(args)){
    factors = readRDS(args[["factors"]])
    if(!"rep_factors" %in% names(args)){
      cat("Using first factor for replicates\n")
      rep_factors = colnames(factors)[1]
    }else{
      rep_factors = args[["rep_factors"]]
    }
  }else{
    factors = NULL
    rep_factors = NULL
  }
  dir.create(path = paste0(wd,"/",Experiment_name),showWarnings = FALSE)
  setwd(paste0(wd,"/",Experiment_name))
  if(sum(c("rds","tsv","csv") %in% names(args)) > 1){
    stop("Error: Select only one data format (rds / tsv / csv)",call. = TRUE)
  }
  if("rds" %in% names(args)){
    cat("Loading ",args[["rds"]],"\n")
    data = readRDS(args[["rds"]])
  }
  if("tsv" %in% names(args)){
    cat("Loading ",args[["tsv"]],"\n")
    data = read.delim(file = args[["tsv"]],header = TRUE,sep = "\t",row.names = 1)
  }
  if("csv" %in% names(args)){
    cat("Loading ",args[["csv"]],"\n")
    data = read.delim(file = args[["csv"]],header = TRUE,sep = ",",row.names = 1)
  }
  threads = as.numeric(args[["threads"]])
  
  ##### Additional arguments #####
  invisible(
    sapply(
      which(names(args) %in% "arg"),
      function(x){
        eval(parse(text = args[[x]]),envir = .GlobalEnv)
      }
    )
  )
  
  cat("Working directory: ",getwd(),"\n",sep = "")
  cat("Experiment name: ",Experiment_name,"\n",sep = "")
  cat("Number of threads: ",threads,"\n",sep = "")

  ##### Register threads #####
  # if(.Platform$OS.type == "unix"){
  #   register(BPPARAM = MulticoreParam(workers = threads))
  # }else{
  #   register(BPPARAM = SerialParam())
  # }
  
  ##### RData output ######
  .classes = NULL
  for(.obj in ls()){
    suppressWarnings({.classes[.obj] = class(get(.obj))})
  }
  prmatrix(matrix(data = c(ls(),.classes),nrow = length(ls()),ncol = 2),quote = FALSE,rowlab = rep("",length(ls())),collab = rep("",2))
  rm(.classes,.obj)
  
}else{
  threads = detectCores()
  wd = rchoose.dir(caption = "Choose working directory:")
  section = rselect.list(choices = c("TOM","power"),multiple = TRUE)
  if("TOM" %in% section){
    softPower = as.numeric(readline(prompt = "Select soft power value"))
    blockwise = rselect.list(choices = c("Single TOM","Blockwise module creation"),multiple = FALSE)
    if(blockwise == "Blockwise module creation"){
      maxBlockSize = readline(prompt = "Select maximum block size: ")
    }
  }
  minModuleSize = readline(prompt = "Select minimum module size (Typically 50): ")
  setwd(wd)
  data_format = rselect.list(choices = c("rds","tsv","csv"),multiple = FALSE,title = "Select data format")
  min_expression = as.numeric(readline(prompt = "Minimum expression: "))
  if(data_format == "rds"){
    data = readRDS(file.choose())
  }
  if(data_format == "tsv"){
    data = read.delim(file = file.choose(),header = TRUE,sep = "\t")
  }
  if(data_format == "csv"){
    data = read.delim(file = file.choose(),header = TRUE,sep = ",")
  }
  factors = readRDS(file.choose())
  rep_factors = select.list(choices = colnames(factors),multiple = TRUE,title = "Select distinct factors for replicates")
}

enableWGCNAThreads(nThreads = threads)

data = t(data[apply(data,1,max) > min_expression,])
if(!is.null(rep_factors)){
  experiment_reps = factor(paste(factors[,rep_factors[1]],sep = "."))
  if(length(rep_factors) > 1){
  for(i in 2:length(rep_factors)){
    experiment_reps = factor(paste(experiment_reps,factors[,rep_factors[i]],sep = "."))
  }
  }
}

if(goodSamplesGenes(data, verbose = 0)$allOK){
  cat("All samples and genes passed QC !\n")
}else{
  cat(paste0("Problematic genes: ",paste(colnames(data)[!goodSamplesGenes(data, verbose = 0)$goodGenes], collapse = ", "),"\n"))
  cat(paste0("Problematic samples: ",paste(rownames(data)[!goodSamplesGenes(data, verbose = 0)$goodSamples],collapse = ", "),"\n"))
  stop("Problem with the data (goodSamplesGenes failed)!", call. = TRUE)
}

if("power" %in% section){
  powers = c(c(1:10), seq(from = 12, to=40, by=2))
  sft = pickSoftThreshold(data,
                        powerVector = powers,
                        verbose = 5,
                        networkType = "signed")

#### Scale Free Topology Model Fit ####
sizeGrWindow(9, 5); par(mfrow = c(1,2)); cex1 = 0.9
pdf(file = paste0(Experiment_name,"_sft.pdf"), width = 9, height = 6)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab=expression("Scale Free Topology Model Fit (signed " ~ R^{2} ~ ")"),type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
while(!is.null(dev.list())) dev.off()
save.image(paste0(Experiment_name,"_sft.RData"))
}


if("TOM" %in% names(args)){
  if(!blockwise){
    dissTOM = 1 - TOMsimilarity(adjacency(data,
                                          power = softPower,
                                          type = "signed"),
                                TOMType="signed")
    
    # save.image(paste0(Experiment_name,"_dissTOM.RData"))
    
    # Tree clustering
    geneTree = flashClust(as.dist(dissTOM), method = "average")
    sizeGrWindow(12,9)
    pdf(file = paste0("Gene_Dendrogram_",Experiment_name,".pdf"), width = 9, height = 6)
    plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
         labels = FALSE, hang = 0.04)
    while(!is.null(dev.list())) dev.off()
    
    # Module construction
    dynamicMods = cutreeDynamic(dendro = geneTree,
                                distM = dissTOM,
                                deepSplit = 2,
                                pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)
    write.table(table(dynamicMods),row.names = FALSE,quote = FALSE, file = paste0("Dynamic_Modules_table_",Experiment_name,".txt"))
    
    # Plotting modules
    dynamicColors = labels2colors(dynamicMods)
    write.table(table(dynamicColors), file = paste0("Modules_table_",Experiment_name,".txt"))
    
    sizeGrWindow(8,6)
    pdf(file = paste0("Gene_dendrogram_modules_",Experiment_name,".pdf"), width = 9, height = 6)
    plotDendroAndColors(dendro = geneTree,
                        colors = dynamicColors,
                        groupLabels = "Dynamic Tree Cut",
                        dendroLabels = FALSE,
                        hang = 0.03,
                        addGuide = TRUE,
                        guideHang = 0.05,
                        main = "Gene dendrogram and module colours")
    while(!is.null(dev.list())) dev.off()
    
    # Calculate modules eigenvalues
    MEList = moduleEigengenes(data, colors = dynamicColors)
    MEDiss = 1-cor(MEList$eigengenes)
    write.csv(MEDiss, file = paste0("MEs_modules_",Experiment_name,".csv"))
    write.csv(MEList$eigengenes, file = paste0("MEs_conditions_",Experiment_name,".csv"))
    
    METree = hclust(as.dist(MEDiss), method = "average")
    sizeGrWindow(7, 6)
    pdf(file = paste0("Dendrogram_module_eigenvalues_",Experiment_name,".pdf"), width = 9, height = 6)
    plot(METree, main = "Clustering module eigenvalues",
         xlab = "", sub = "")
    while(!is.null(dev.list())) dev.off()
    
  }else{
    blockTOM = blockwiseModules(data,
                                maxBlockSize = maxBlockSize,
                                power = softPower,
                                TOMType = "signed",
                                minModuleSize = minModuleSize,
                                reassignThreshold = 0,
                                mergeCutHeight = 0.25,
                                numericLabels = TRUE,
                                saveTOMs = TRUE,
                                saveTOMFileBase = Experiment_name,
                                verbose = 3)
    
    # save.image(paste0(Experiment_name,"_dissTOM.RData"))
    module_colors = labels2colors(blockTOM$colors)
    names(module_colors) = names(blockTOM[["colors"]])
    
    for(i in 1:length(blockTOM$dendrograms)){
      png(filename = paste0("Module_tree_",Experiment_name,"_",i,".png"), width = 1920, height = 1080, units = "px")
      plotDendroAndColors(blockTOM$dendrograms[[i]],
                          module_colors[blockTOM$blockGenes[[i]]],
                          "Module colors",
                          main = paste0("Gene dendrogram and module colors in block ",i),
                          dendroLabels = FALSE,
                          hang = 0.03,
                          addGuide = TRUE,
                          guideHang = 0.05)
      while(!is.null(dev.list())) dev.off()
    }
    
    MEList = moduleEigengenes(data, colors = module_colors)
    METree = hclust(as.dist(1-cor(MEList$eigengenes)), method = "average")
    png(filename = paste0("Eigengene_module_clustering",Experiment_name,".png") ,width = 1920, height = 1080, units = "px")
    plot(METree,
         main = "Clustering module eigenvalues",
         xlab = "",
         sub = "")
    while(!is.null(dev.list())) dev.off()
  }
  save(list = ls()[grep(pattern = "dissTOM|blockTOM",x = ls(),invert = TRUE)],file = paste0(Experiment_name,"_TOM.RData"))
  # save.image(paste0(Experiment_name,"_TOM.RData"))

ME_rep = function(ME_data){
  return(list(Mean = sapply(levels(experiment_reps),function(x){
    mean(ME_data[experiment_reps %in% x])
  }),
  SD = sapply(levels(experiment_reps),function(x){
    sd(ME_data[experiment_reps %in% x])
  })))
}

ME_eigen_rep = lapply(MEList$eigengenes,ME_rep)
ME_averageExpr_rep = lapply(MEList$averageExpr,ME_rep)

# Create eigen plots
dir.create("eigen_plots",showWarnings = FALSE)
for(me in names(MEList$eigengenes)){
  png(filename = paste0("eigen_plots/",me,".png"),width = 1920,height = 1080,units = "px")
  print(ggplot(data = data.frame(Sample = names(ME_eigen_rep[[me]][["Mean"]]),Mean = ME_eigen_rep[[me]][["Mean"]]),aes(x = Sample,y = Mean,group=1)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin=ME_eigen_rep[[me]][["Mean"]]-ME_eigen_rep[[me]][["SD"]],
                    ymax=ME_eigen_rep[[me]][["Mean"]]+ME_eigen_rep[[me]][["SD"]]),
                width=.2,
                position=position_dodge(0.05)))
  while (!is.null(dev.list())){dev.off()}
}

# Create Average expression plots
dir.create("Average_Expr_plots",showWarnings = FALSE)
for(me in names(MEList$averageExpr)){
  png(filename = paste0("Average_Expr_plots/",me,".png"),width = 1920,height = 1080,units = "px")
  print(ggplot(data = data.frame(Sample = names(ME_averageExpr_rep[[me]][["Mean"]]),Mean = ME_averageExpr_rep[[me]][["Mean"]]),aes(x = Sample,y = Mean,group=1)) +
          geom_line() +
          geom_point() +
          geom_errorbar(aes(ymin=ME_averageExpr_rep[[me]][["Mean"]]-ME_averageExpr_rep[[me]][["SD"]],
                            ymax=ME_averageExpr_rep[[me]][["Mean"]]+ME_averageExpr_rep[[me]][["SD"]]),
                        width=.2,
                        position=position_dodge(0.05)))
  while (!is.null(dev.list())){dev.off()}
}
}

# Top membership genes
geneModuleMembership = as.data.frame(cor(data, orderMEs(moduleEigengenes(data, dynamicColors)$eigengenes), use = "p"))
topgenes = lapply(colnames(geneModuleMembership),function(module){
  temp = geneModuleMembership[[module]]
  names(temp) = rownames(geneModuleMembership)
  return(sort(temp,decreasing = TRUE))
})
names(topgenes) = colnames(geneModuleMembership)
dir.create("module_genes",showWarnings = FALSE)
for(module in names(topgenes)){
  write.table(as.data.frame(topgenes[[module]]),
              file = paste0("module_genes/",module,".txt"),
              quote = FALSE,
              sep = "\t",
              row.names = TRUE,
              col.names = FALSE)
}


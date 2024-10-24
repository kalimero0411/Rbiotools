##### SETUP #####
packages = c("purrr", "bindr", "Biobase", "BiocGenerics", "checkmate", "digest","dplyr", "dynamicTreeCut", "flashClust", 
             "ggplot2", "gplots", "ggpubr", "IRanges", "RColorBrewer", "robust", "tidyr", "WGCNA","parallel","R.utils",
             "MatrixGenerics","pheatmap")

if(interactive()){
  packages = c(packages,"BiocParallel")
}

invisible(
  suppressMessages(
    if(!require("BiocManager",character.only = TRUE,quietly = TRUE)){
      cat("Installing BiocManager\n",sep = "")
      install.packages("BiocManager")
    }))

cat("#####   Loading packages   #####\n")
invisible(
  suppressMessages(
    lapply(packages,function(x){
      if(!require(x,character.only = TRUE,quietly = TRUE)){
        cat("Installing package: ",x,"\n",sep = "")
        BiocManager::install(x,update = FALSE,ask = FALSE)
        library(x,character.only = TRUE,quietly = TRUE)
      }
    })))


options(stringsAsFactors = FALSE)

if(!interactive()){
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  must_args = c("rdata","wd","name","threads")
  if(!all(must_args %in% names(args))){
    help = matrix(data = c("--rdata    ","RData file path",
                           "--wd","Working directory path",
                           "--name","Experiment name",
                           "--factors","Experimantal factors to test(Seperated by comma; Default = All factors)",
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
    stop(paste(must_args[!must_args %in% names(args)],collapse = " | "), call. = TRUE)
  }
  
  # saveRDS(data,file = "data.rds")
  cat("Loading RData file: ",args[["rdata"]],"\n",sep = "")
  if("rdata" %in% names(args)){
    load(args[["rdata"]])
  }
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  section = c()
  wd = args[["wd"]]
  Experiment_name = args[["name"]]
  if("factors" %in% names(args)){
    rep_factors = unlist(strsplit(args[["factors"]],split = ","))
  }else{
    rep_factors = colnames(experimental_design)
  }
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
  dir.create(path = paste0(wd,"/",Experiment_name),showWarnings = FALSE)
  setwd(paste0(wd,"/",Experiment_name))
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
  wd = rstudioapi::selectDirectory(caption = "Choose working directory:")
  load(file.choose())
  Experiment_name = as.character(readline(prompt = "Select experiment name: "))
  section = select.list(choices = c("TOM","power"),multiple = TRUE)
  if("TOM" %in% section){
    softPower = as.numeric(readline(prompt = "Select soft power value: "))
    blockwise = "Blockwise module creation" == select.list(choices = c("Single TOM","Blockwise module creation"),multiple = FALSE)
    if(blockwise){
      maxBlockSize = readline(prompt = "Select maximum block size: ")
    }
  }
  minModuleSize = as.numeric(readline(prompt = "Select minimum module size (Typically 50): "))
  setwd(wd)
  data_format = select.list(choices = c("rds","tsv","csv"),multiple = FALSE,title = "Select data format")
  min_expression = as.numeric(readline(prompt = "Minimum expression: "))
  rep_factors = select.list(choices = colnames(experimental_design),multiple = TRUE,title = "Select factors")
}

enableWGCNAThreads(nThreads = threads)

if(!exists("input_mat")){
  input_mat = t(data_vst[apply(data_vst,1,max) > min_expression,])
if(!is.null(rep_factors)){
  experiment_reps = factor(paste(experimental_design[,rep_factors[1]],sep = "."))
  if(length(rep_factors) > 1){
  for(i in 2:length(rep_factors)){
    experiment_reps = factor(paste(experiment_reps,experimental_design[,rep_factors[i]],sep = "."))
  }
  }
}
}

if(goodSamplesGenes(input_mat, verbose = 0)$allOK){
  cat("All samples and genes passed QC !\n")
}else{
  cat(paste0("Problematic genes: ",paste(colnames(input_mat)[!goodSamplesGenes(input_mat, verbose = 0)$goodGenes], collapse = ", "),"\n"))
  cat(paste0("Problematic samples: ",paste(rownames(input_mat)[!goodSamplesGenes(input_mat, verbose = 0)$goodSamples],collapse = ", "),"\n"))
  stop("Problem with the data (goodSamplesGenes failed)!", call. = TRUE)
}

if("power" %in% section){
  powers = c(c(1:10), seq(from = 12, to=40, by=2))
  sft = pickSoftThreshold(input_mat,
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
abline(h=0.80,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
while(!is.null(dev.list())) dev.off()
save.image(paste0(Experiment_name,"_sft.RData"))
}


if("TOM" %in% section){
  if(!blockwise){
    dissTOM = 1 - TOMsimilarity(adjacency(input_mat,
                                          power = softPower,
                                          type = "signed"),
                                TOMType="signed")
    
    save.image(paste0(Experiment_name,"_dissTOM.RData"))
    
    # Tree clustering
    geneTree = flashClust(as.dist(dissTOM), method = "average")
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
    names(dynamicColors) = colnames(input_mat)
    write.table(file = paste0("Gene_module_map_",Experiment_name,".txt"),data.frame(Gene = colnames(input_mat),Module = dynamicColors),col.names = TRUE,quote = FALSE,row.names = FALSE,sep = "\t")
    write.table(table(dynamicColors), file = paste0("Modules_table_",Experiment_name,".txt"))
    
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
    MEList = moduleEigengenes(input_mat, colors = dynamicColors)
    MEDiss = 1-cor(MEList$eigengenes)
    write.csv(MEDiss, file = paste0("MEs_modules_",Experiment_name,".csv"))
    write.csv(MEList$eigengenes, file = paste0("MEs_conditions_",Experiment_name,".csv"))
    
    METree = hclust(as.dist(MEDiss), method = "average")
    pdf(file = paste0("Dendrogram_module_eigenvalues_",Experiment_name,".pdf"), width = 9, height = 6)
    plot(METree, main = "Clustering module eigenvalues",
         xlab = "", sub = "")
    while(!is.null(dev.list())) dev.off()
    
  }else{
    blockTOM = blockwiseModules(input_mat,
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
    
    MEList = moduleEigengenes(input_mat, colors = module_colors)
    METree = hclust(as.dist(1-cor(MEList$eigengenes)), method = "average")
    png(filename = paste0("Eigengene_module_clustering",Experiment_name,".png") ,width = 1920, height = 1080, units = "px")
    plot(METree,
         main = "Clustering module eigenvalues",
         xlab = "",
         sub = "")
    while(!is.null(dev.list())) dev.off()
  }


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
save(list = ls()[grep(pattern = "dissTOM|blockTOM",x = ls(),invert = TRUE)],file = paste0(Experiment_name,"_TOM.RData"))
# save.image(paste0(Experiment_name,"_TOM.RData"))
}

if("top genes" %in% section){
# Top membership genes
geneModuleMembership = as.data.frame(cor(input_mat, orderMEs(moduleEigengenes(input_mat, dynamicColors)$eigengenes), use = "p"))
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

t.test2 = function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}

Aemods = c("blue",
       "brown",
       "darkgreen",
       "darkgrey",
       "darkturquoise",
       "lightgreen",
       "purple",
       "salmon",
       "turquoise")

Aemods = c("darkmagenta",
           "darkgrey",
           "cyan",
           "royalblue",
           "salmon",
           "tan",
           "violet",
           "yellow",
           "black")

top_stats = lapply(Aemods,function(Ae){
  res = lapply(names(topgenes[[paste0("ME",Ae)]][1:100]),function(topgene){
  top = t(TMM_non_zero[topgene,])[,1]
  top_stats = sapply(levels(experimental_design[[1]]),function(fac){
    samples = rownames(experimental_design)[experimental_design[[1]] %in% fac]
    return(c(
      Mean = mean(top[names(top) %in% samples]),
      SD = sd(top[names(top) %in% samples]),
      N = length(top[names(top) %in% samples])))
  })
  return(top_stats)
  })
  names(res) = names(topgenes[[paste0("ME",Ae)]][1:100])
  return(res)
})
names(top_stats) = Aemods

top_tests = lapply(Aemods,function(Ae){
tests = lapply(top_stats[[Ae]],function(stats){
  sapply(paste0("T",c(0,1,2,5,8)),function(x){
t.test2(m1 = stats["Mean",paste0(x,".31")],
        m2 = stats["Mean",paste0(x,".27")],
        s1 = stats["SD",paste0(x,".31")],
        s2 = stats["SD",paste0(x,".27")],
        n1 = stats["N",paste0(x,".31")],
        n2 = stats["N",paste0(x,".27")],
        equal.variance = FALSE)[["p-value"]]
})
})
})

names(top_tests) = Aemods

sig_pattern = unname(as.data.frame(t(expand.grid(0,c(0,1),c(0,1),1,c(0,1)))))

top_sig_pattern = lapply(sig_pattern,function(pattern){
  top = lapply(top_tests,function(x){
  names(x)[sapply(x,function(y){
    return(identical(unname(y)<0.05,as.logical(pattern)))
  })]
})
  temp = unlist(top,use.names = FALSE)
  temp = setNames(temp,rep(names(top),lengths(top)))
  temp = temp[!duplicated(temp)]
  return(temp)
  })
top_sig_pattern = unlist(top_sig_pattern)
top_sig_stats = apply(t(TMM_non_zero[top_sig_pattern,]),MARGIN = 2,ME_rep)

# experimental_design = data.frame(Timepoint = rep(paste0("T",c(0,1,2,5,8)),each = 2),
#                                  Temperature = rep(c(27,31),times = 5))
dir.create("module_genes/top_genes",showWarnings = FALSE)
invisible(sapply(names(top_sig_stats),function(me){
  top_data = top_sig_stats[[me]]
  top_control = data.frame(Sample = unique(gsub(names(top_data[["Mean"]]),pattern = "[.].+",replacement = "",perl = TRUE)),
                           Mean = top_data[["Mean"]][grepl(pattern = "27",x = names(top_data[["Mean"]]))],
                           SD = top_data[["SD"]][grepl(pattern = "27",x = names(top_data[["SD"]]))])
  top_test = data.frame(Sample = unique(gsub(names(top_data[["Mean"]]),pattern = "[.].+",replacement = "",perl = TRUE)),
                           Mean = top_data[["Mean"]][grepl(pattern = "31",x = names(top_data[["Mean"]]))],
                           SD = top_data[["SD"]][grepl(pattern = "31",x = names(top_data[["SD"]]))])
  png(filename = paste0("module_genes/top_genes/",me,".png"),width = 1080,height = 1080,units = "px")
  print(ggplot(data = top_control, aes(x = Sample,y = Mean,group = 1)) +
          geom_errorbar(aes(ymin = Mean - SD,
                            ymax = Mean + SD),
                        width=.2,
                        position=position_dodge(0.05),
                        color = "black") +
          geom_line(aes(y = Mean)) +
          geom_point(aes(y = Mean), color = "black") +
          geom_line(data = top_test, aes(x = Sample,y = Mean,group = 1),color = "red") +
          geom_errorbar(data = top_test,
                        aes(ymin = Mean - SD,
                            ymax = Mean + SD),
                        width=.2,
                        position=position_dodge(0.05),
                        color = "red") +
          geom_point(data = top_test,aes(y = Mean), color = "red") +
          ggtitle(paste0(me," (",names(top_sig_pattern[top_sig_pattern %in% me]),")")) +
          ylab("Mean TMM") +
          theme(text = element_text(size = 30)))
  while (!is.null(dev.list())){dev.off()}
}))

write.table(names(top_sig_stats),file = "top_genes_ttest.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

  data_topgenes_scaled = scale(data_vst)[rownames(data_vst) %in% names(top_sig_stats),]

data_topgenes_scaled = data_topgenes_scaled[,order(match(coldata[[1]],levels(coldata[[1]]))),drop = FALSE]
coldata = coldata[order(match(coldata[[1]],levels(coldata[[1]]))),,drop = FALSE]

pheatmap(data_topgenes_scaled,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         annotation_col = coldata,
         filename = "top_genes.png")

}

#### DEG heatmaps

AE = t(moduleEigengenes(input_mat[,colnames(input_mat) %in% degs_all,drop = FALSE], colors = dynamicColors[names(dynamicColors) %in% degs_all])[["averageExpr"]])
rownames(AE) = gsub(pattern = "^AE",replacement = "",x = rownames(AE))
png(filename = paste0("DEG_heamap_all_",Experiment_name,".png"), width = 720, height = 1080, units = "px")
pheatmap(AE,annotation_col = experimental_design,main = "DEG modules (all changes)",fontsize = 20)
while(!is.null(dev.list())) dev.off()

AE = t(moduleEigengenes(input_mat[,colnames(input_mat) %in% degs_up,drop = FALSE], colors = dynamicColors[names(dynamicColors) %in% degs_up])[["averageExpr"]])
rownames(AE) = gsub(pattern = "^AE",replacement = "",x = rownames(AE))
png(filename = paste0("DEG_heamap_up_",Experiment_name,".png"), width = 720, height = 1080, units = "px")
pheatmap(AE,annotation_col = experimental_design,main = "DEG modules (upregulated)",fontsize = 20)
while(!is.null(dev.list())) dev.off()

AE = t(moduleEigengenes(input_mat[,colnames(input_mat) %in% degs_down,drop = FALSE], colors = dynamicColors[names(dynamicColors) %in% degs_down])[["averageExpr"]])
rownames(AE) = gsub(pattern = "^AE",replacement = "",x = rownames(AE))
png(filename = paste0("DEG_heamap_down_",Experiment_name,".png"), width = 720, height = 1080, units = "px")
pheatmap(AE,annotation_col = experimental_design,main = "DEG modules (downregulated)",fontsize = 20)
while(!is.null(dev.list())) dev.off()


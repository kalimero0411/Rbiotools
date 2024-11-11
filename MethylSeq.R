#######     DNA methylation analysis      ##########

packages=c("methylKit","Rsamtools","genomation","goseq","BiocParallel","parallel","tools","factoextra",
           "R.utils","ggplot2","rstudioapi","matrixStats")

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
init_params = list()

###### Cluster commands ######
if(!interactive()){
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  must_args = c("wd","name")
  if(!all(must_args %in% names(args))){
    help = matrix(data = c("--wd","Working directory",
                    "--name","Experiment name (all data will output to a directory by that name in the working directory)",
                    "--bismark","Bismark output BAM files",
                    "--analyze","Skip reading input and Analyze data from RData file",
                    "--read","Only read data and output RData file (don't analyze)",
                    "--anno","Annotation file (12 column bed format; convert GFF/GTF with agat_convert_sp_gff2bed.pl)",
                    "--exp","Experimental design file",
                    "--control","Control variable for each factor (comma separated; in the same order as in --exp; Default = first variable in each factor)",
                    "--context","Methylation context: CpG, CHG and CHH (Default = All)",
                    "--t","Number of compute threads",
                    "--arg","Additional R arguments (multiple arguments in separate flags)")
                  ,ncol = 2,byrow = TRUE)
    prmatrix(help,quote = FALSE,rowlab = rep("",nrow(help)),collab = rep("",2))
    stop(paste0("Missing command line input --> ",paste(must_args[!must_args %in% names(args)],collapse = " | ")), call. = TRUE)
  }
  
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  init_params[["wd"]] = normalizePath(args[["wd"]])
  init_params[["name"]] = args[["name"]]
  
  if("read" %in% names(args) == "analyze" %in% names(args)){
    init_params[["section"]] = c("Read","Analyze")
  }else{
    init_params[["section"]] = c("Read","Analyze")[c("read" %in% names(args),"analyze" %in% names(args))]
    if("Analyze" %in% init_params[["section"]]){
      init_params[["analyze"]] = normalizePath(args[["analyze"]])
    }
  }
  if("Read" %in% init_params[["section"]]){
    init_params[["exp"]] = normalizePath(args[["exp"]])
    experimental_design = read.table(file = init_params[["exp"]],header = TRUE,sep = "\t")
    experimental_design[["File"]] = normalizePath(experimental_design[["File"]])
    if("control" %in% names(args)){
      init_params[["controls"]] = strsplit(x = gsub(pattern = " ",replacement = "",args[["control"]]),split = ",")[[1]]
    }else{
      init_params[["controls"]] = experimental_design[1,3:ncol(experimental_design)]
    }
    names(init_params[["controls"]]) = colnames(experimental_design)[3:ncol(experimental_design)]
    init_params[["annotation"]] = normalizePath(args[["anno"]])
  }
  
  if("context" %in% names(args)){
    init_params[["context"]] = c("CpG","CHG","CHH")[sapply(c("cpg","chg","chh"),function(x){
      any(grepl(pattern = paste0("^",x),strsplit(gsub(pattern = " ",replacement = "",args[["context"]]),split = ",")[[1]],ignore.case = TRUE))
    })]
  }else{
    init_params[["context"]] = c("CpG","CHG","CHH")
  }
  
  if("t" %in% names(args)){
    init_params[["threads"]] = as.numeric(args[["t"]])
  }else{
    init_params[["threads"]] = detectCores()
  }
    
  cat("Number of threads: ",init_params[["threads"]],"\n", sep = "")
  
  ##### Additional arguments #####
  if("arg" %in% names(args)){
  add_args = unname(args[which(names(args) %in% "arg")])
  for(i in seq_along(add_args)){
    cat("Performing: ",add_args[[i]],"\n", sep = "")
    eval(parse(text = add_args[[i]]),envir = .GlobalEnv)
  }
  }
  
##### RData output ######
.classes = NULL
for(.obj in ls()){
  suppressWarnings({.classes[.obj] = class(get(.obj))})
}
prmatrix(matrix(data = c(ls(),.classes),nrow = length(ls()),ncol = 2),quote = FALSE,rowlab = rep("",length(ls())),collab = rep("",2))
rm(.classes,.obj)

}else{
  
  init_params[["threads"]] = detectCores()
  
  ###### Set working directory ######
  init_params[["wd"]] = selectDirectory(caption = "Select working directory")
  init_params[["name"]] = readline(prompt = "Experiment name: ")
  init_params[["section"]] = select.list(choices = c("Read","Analyze"),multiple = TRUE,title = "Skip Reading data?",preselect = c("Read","Analyze"))
  if("Analyze" %in% init_params[["section"]] & !"Read" %in% init_params[["section"]]){
    init_params[["analyze"]] = selectFile(caption = "Select RData file to analyze")
  }else{
    ###### Input analysis settings ######
    cat("#####   Input analysis settings   #####\n")

    init_params[["exp"]] = normalizePath(selectFile(caption = "Select experimental design file",path = init_params[["wd"]]))
    experimental_design = read.table(file = init_params[["exp"]],header = TRUE,sep = "\t")
    experimental_design[["File"]] = normalizePath(experimental_design[["File"]])
    init_params[["controls"]] = c()
    for(i in 3:ncol(experimental_design)){
      init_params[["controls"]] = c(init_params[["controls"]],select.list(unique(experimental_design[[i]]),multiple = FALSE,title = "Select control variable",graphics = TRUE))
    }
    names(init_params[["controls"]]) = colnames(experimental_design[3:ncol(experimental_design)])
    
    init_params[["annotation"]] = normalizePath(selectFile(caption = "Select genome annotation file",path = init_params[["wd"]]))
    init_params[["assembly"]] = readline(prompt = "Organism: ")
    init_params[["context"]] = select.list(choices = c("CpG","CHG","CHH"),multiple = TRUE,title = "Select contexts to analyze",graphics = TRUE,preselect = c("CpG","CHG","CHH"))
    
    init_params[["window_size"]] = as.numeric(readline(prompt = "Window size for differentially methylated regions (default = 100): "))
    init_params[["alpha"]] = as.numeric(readline(prompt = "pvalue cutoff for methylation difference (default = 0.05): "))
    init_params[["methylation_diff"]] = as.numeric(readline(prompt = "Methylation difference cutoff percent (default = 25): "))
    if(is.na(init_params[["window_size"]])){init_params[["window_size"]] = 100}
    if(is.na(init_params[["alpha"]])){init_params[["alpha"]] = 0.05}
    if(is.na(init_params[["methylation_diff"]])){init_params[["methylation_diff"]] = 25}
    
  }
}

start_time=Sys.time()

cat("Working directory: ",getwd(),"\n", sep = "")
cat("Experiment name: ",init_params[["name"]],"\n", sep = "")
setwd(init_params[["wd"]])

if("analyze" %in% names(init_params)){
  init_params_rem = init_params
  load(init_params[["analyze"]])
  for(i in names(init_params_rem)){
    init_params[[i]] = init_params_rem[[i]]
  }
  rm(init_params_rem)
}

##### Register threads #####
if(.Platform$OS.type == "unix"){
  register(BPPARAM = MulticoreParam(workers = init_params[["threads"]]))
}else{
  init_params[["threads"]] = 1
  register(BPPARAM = SerialParam())
}

if(!exists("genome_annotation")){
  genome_annotation = readTranscriptFeatures(location = init_params[["annotation"]])
}

for(i in names(init_params[["controls"]])){
  init_params[["factors"]][[i]] = as.numeric(relevel(factor(experimental_design[[i]]),ref = init_params[["controls"]][i])) - 1
  names(init_params[["factors"]][[i]]) = relevel(factor(experimental_design[[i]]),ref = init_params[["controls"]][i])
}

#####   Run files in selected contexts   ######

#####  Import samples  #####
if("Read" %in% init_params[["section"]]){
  dir.create("Basic_stats",showWarnings = FALSE)
  meth_unite = list()
  meth_unite_feature = list()
  
  for(variable in names(init_params[["factors"]])){
  dir.create(paste0("MethylDB/",variable),showWarnings = FALSE,recursive = TRUE)
  processBismarkAln(location = as.list(experimental_design[["File"]]),
                    sample.id = as.list(experimental_design[["Sample"]]),
                    assembly = init_params[["assembly"]],
                    treatment = init_params[["factors"]][[variable]],
                    save.context = init_params[["context"]],
                    save.folder = paste0("MethylDB/",variable))
  for(context in init_params[["context"]]){
    init_params[["MethylDB"]][[variable]][[context]] = normalizePath(path = paste0("MethylDB/",variable,"/",experimental_design[["Sample"]],"_",context,".txt"))
  }
  }
  
  for(variable in names(init_params[["factors"]])){
  for(context in init_params[["context"]]){
  cat("#####   Importing sample files in ",context," context   #####\n", sep = "")
  raw_meth = methRead(location = as.list(init_params[["MethylDB"]][[variable]][[context]]),
                      sample.id = as.list(experimental_design[["Sample"]]),
                      pipeline = "bismark",
                      treatment = init_params[["factors"]][[variable]],
                      assembly = init_params[["assembly"]],
                      context = context)
  names(raw_meth) = experimental_design[["Sample"]]

      ######  Basic stats in numbers  ######
    cat("#####   Calculating basic stats, plots and coverage plots   #####\n")
      for(i in experimental_design[["Sample"]]) {
        methylstats[[variable]][[context]][[i]] = capture.output(getMethylationStats(raw_meth[[i]],plot = FALSE,both.strands = FALSE),file = paste0("Basic_stats/Methyl_stats_",i,"_",variable,"_",context,".methylstats"))
        
        ###### Basic stats plots ######
        png(filename = paste0("Basic_stats/Methyl_stats_",i,"_",variable,"_",context,".png"),width = 1920,height = 1080,units = "px")
        getMethylationStats(raw_meth[[i]],plot = TRUE,both.strands = FALSE)
        while (!is.null(dev.list())) dev.off()
        
        ###### Coverage stats ######
        png(filename = paste0("Basic_stats/Coverage_stats_",i,"_",variable,"_",context,".png"),width = 1920,height = 1080,units = "px")
        getCoverageStats(raw_meth[[i]],plot = TRUE,both.strands = FALSE)
        while (!is.null(dev.list())) dev.off()
      }
      
      ###### Filter and normalize each sample by read coverage ######
      cat("#####   Filtering data   #####\n")
      raw_meth_filter = filterByCoverage(methylObj = raw_meth,
                                         lo.count = 10,
                                         lo.perc = NULL,
                                         hi.count = NULL,
                                         hi.perc = 99.9
                                         )
      rm(raw_meth)
      gc()
      cat("#####   Normalizing data   #####\n")
      raw_meth_norm = normalizeCoverage(obj = raw_meth_filter,
                                        method = "median"
                                        )
      
      rm(raw_meth_filter)
      gc()
      ###### Create tiles ######
      cat("#####   Creating region tiles   #####\n")
      tiles = tileMethylCounts(raw_meth_norm,
                               win.size = init_params[["window_size"]],
                               step.size = init_params[["window_size"]],
                               mc.cores = init_params[["threads"]]
                               )
      
      ###### Unite samples ######
      for(scope in c("SMP","DMR")){
      cat("#####   Uniting samples   #####\n")
        meth_unite[[scope]][[variable]][[context]] = unite(object = if(scope == "SMP"){raw_meth_norm}else{tiles},destrand = context == "CpG",save.db = FALSE,mc.cores = init_params[["threads"]])
      
      #### Filter low SD ####
        cat("#####   Filtering low SD sites    #####\n")
        meth_unite[[scope]][[variable]][[context]] = meth_unite[[scope]][[variable]][[context]][rowSds(percMethylation(meth_unite[[scope]][[variable]][[context]]),na.rm = TRUE) > 2]
        
      ###### Data lists ######

        ###### Extract regions ######
        for(feature in c("promoters","exons","introns","TSSes")){
          cat("#####   Extracting regions for ",feature,"   #####\n", sep = "")
          meth_unite_feature[[scope]][[variable]][[context]][[feature]] = regionCounts(meth_unite[[scope]][[variable]][[context]],
                                                            regions = genome_annotation[[feature]],
                                                            mc.cores = init_params[["threads"]])
        }
      }
      rm(raw_meth_norm,tiles)
      gc()
      
  }
  }
  save.image(paste0(init_params[["name"]],"_imported_data.RData"))
  }
  
  if("Analyze" %in% init_params[["section"]]){
  cat("#####   Including following contexts in analysis: ",paste0(init_params[["context"]],collapse = ", "),"   #####\n", sep = "")

    perc_meth = list()
    meth_poly = list()
    meth_diff = list()
    meth_genes = list()
    meth_distTSS = list()
    
    dir.create("Raw_data",showWarnings = FALSE)
    dir.create("Correlation",showWarnings = FALSE)
    dir.create("Cluster",showWarnings = FALSE)
    dir.create("PCA",showWarnings = FALSE)
    dir.create("Diff_meth",showWarnings = FALSE)
    dir.create("Annotation",showWarnings = FALSE)
    dir.create("Volcano",showWarnings = FALSE)
    
  ###### Annotation and Regions loop ######
  for(variable in names(init_params[["factors"]])){
  for(context in init_params[["context"]]){
  for(scope in c("SMP","DMR")){
  cat("######  Running analysis for ",variable," in ",context," context for ",scope,"s ######\n", sep = "")
    
    if(nrow(meth_unite[[scope]][[variable]][[context]]) > 0){
      write.table(x = getData(meth_unite[[scope]][[variable]][[context]]),
                  file = paste0("Raw_data/",scope,"_",variable,"_",context,".txt"),
                  quote = FALSE,
                  sep = "\t",
                  row.names = FALSE,
                  col.names = TRUE)
      perc_meth[[scope]][[variable]][[context]] = list(mean = colMeans(percMethylation(methylBase.obj = meth_unite[[scope]][[variable]][[context]])),
                                                       sd = colSds(percMethylation(methylBase.obj = meth_unite[[scope]][[variable]][[context]])))
    }
  
  if(!any(colSums(getData(meth_unite[[scope]][[variable]][[context]])[,grep(pattern = "numCs",meth_unite[[scope]][[variable]][[context]]@names)]) == 0)){
    
  ###### Correlation ######
  cat("#####   Creating correlation stats and plots   #####\n")
  
  png(paste0("Correlation/Correlation_",scope,"_",variable,"_",context,".png"),width = 1920,height = 1080,units = "px")
  capture.output(getCorrelation(meth_unite[[scope]][[variable]][[context]],plot = TRUE),file = paste0("Correlation/Correlation_",scope,"_",variable,"_",context,".txt"))
  while (!is.null(dev.list())) dev.off()
  
  ###### Cluster the samples ######
  cat("#####   Creating cluster plots   #####\n")
  png(filename = paste0("Cluster/Cluster_",scope,"_",variable,"_",context,".png"),width = 1080,height = 1080,units = "px")
  invisible(clusterSamples(meth_unite[[scope]][[variable]][[context]],dist = "correlation",method = "ward.D2",plot = TRUE))
  while (!is.null(dev.list())) dev.off()
  
  ###### PCA variances and plots ######
  cat("#####   Creating PCA variances and plots   #####\n")
  png(filename = paste0("PCA/ScreePlot_",scope,"_",variable,"_",context,".png"),width = 1920,height = 1080,units = "px")
  PCASamples(meth_unite[[scope]][[variable]][[context]], screeplot=TRUE)
  while (!is.null(dev.list())) dev.off()
  
  png(filename = paste0("PCA/PCA_",scope,"_",variable,"_",context,".png"),width = 1920,height = 1080,units = "px")
  print(fviz_pca_ind(X = PCASamples(meth_unite[[scope]][[variable]][[context]],obj.return = TRUE),repel = TRUE,habillage = experimental_design[[variable]],title = paste0("PCA"),labelsize = 8, pointsize = 3,legend.title = variable) +
          theme(title = element_text(size = 20),axis.title = element_text(size = 20),axis.text = element_text(size = 20),legend.text = element_text(size = 20),legend.title = element_text(size = 22)))
  while (!is.null(dev.list())){dev.off()}
  
  for(variable_idx in sort(unique(init_params[["factors"]][[variable]]))[-1]){
    comp_var = c(unique(names(init_params[["factors"]][[variable]])[init_params[["factors"]][[variable]] == 0]),
                 unique(names(init_params[["factors"]][[variable]])[init_params[["factors"]][[variable]] == variable_idx]))
  
  ###### Calculate Single methylation polymorphisms (SMPs) ######
    cat("#####   Subsetting into samples ",comp_var[1]," and ",comp_var[2]," for ",scope,"s  #####\n", sep = "")
    subset_poly = reorganize(methylObj = meth_unite[[scope]][[variable]][[context]],
                            sample.ids = experimental_design[["Sample"]][init_params[["factors"]][[variable]] %in% c(0,variable_idx)],
                            treatment = init_params[["factors"]][[variable]][init_params[["factors"]][[variable]] %in% c(0,variable_idx)])
    
    save.image("tmp6.RData")
    stop("OK stop")

    meth_poly[[scope]][[variable]][[context]][[comp_var[2]]] = calculateDiffMeth(subset_poly,overdispersion = "MN",mc.cores = init_params[["threads"]],save.db = FALSE)

    volcano_data = getData(meth_poly[[scope]][[variable]][[context]][[comp_var[2]]])[,c("meth.diff","qvalue")]
    volcano_data[["Expression"]] = apply(volcano_data,MARGIN = 1,function(x){
      if(abs(as.numeric(x[1])) >= init_params[["methylation_diff"]] & as.numeric(x[2]) < init_params[["alpha"]]){
        if(as.numeric(x[1]) > 0){
          "Hyper-methylated"
        }else{
          "Hypo-methylated"
        }
      }else{
        "Unchanged"
      }
    })

    volcano_plot = ggplot(volcano_data, aes(meth.diff, -log(qvalue,10))) +
      geom_point(aes(color = Expression,alpha = Expression), size = 3) +
      xlab(expression("Methylation difference")) +
      ylab(expression(paste("-log"[10]*"(",italic("p"),"-value)"))) +
      scale_color_manual(values = c("Hyper-methylated" = "dodgerblue3", "Unchanged" = "gray50", "Hypo-methylated" = "firebrick3")) +
      scale_alpha_manual(values = c("Hyper-methylated" = 1,"Unchanged" = 0.25,"Hypo-methylated" = 1)) +
      guides(colour = guide_legend(override.aes = list(size = 3))) +
      geom_hline(yintercept = -log(init_params[["alpha"]],10), linetype="dashed") +
      geom_vline(xintercept = c(-init_params[["methylation_diff"]],init_params[["methylation_diff"]]), linetype="dashed") +
      theme(text = element_text(size = 30))
    
    png(filename = paste0("Volcano/Volcano_",scope,"_",variable,"_",context,"_",comp_var[2],".png"),width = 1920,height = 1080,units = "px")
    print(volcano_plot)
    while (!is.null(dev.list())) dev.off()
    
    rm(volcano_data,volcano_plot)
    
  ###### Get differential methylation stats ######
  cat("#####   Creating stats and plots for SMPs and DMRs for chromosomes   #####\n")
  
  if(sum(meth_poly[[scope]][[variable]][[context]][[comp_var[2]]][["qvalue"]] <= init_params[["alpha"]] & abs(meth_poly[[scope]][[variable]][[context]][[comp_var[2]]][["meth.diff"]]) >= init_params[["methylation_diff"]]) > 0){
    write.table(diffMethPerChr(meth_poly[[scope]][[variable]][[context]][[comp_var[2]]],
                               plot = FALSE,
                               qvalue.cutoff = init_params[["alpha"]],
                               meth.cutoff = init_params[["methylation_diff"]])[["diffMeth.per.chr"]],
                file = paste0("Diff_meth/Diff_meth_per_chr_",scope,"_",variable,"_",context,"_",comp_var[2],".txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
      
    png(filename = paste0("Diff_meth/Diff_meth_per_chr_",scope,"_",variable,"_",context,"_",comp_var[2],".png"),width = 1920,height = 1080,units = "px")
    diffMethPerChr(meth_poly[[scope]][[variable]][[context]][[comp_var[2]]],
                   plot = TRUE,
                   qvalue.cutoff = init_params[["alpha"]],
                   meth.cutoff = init_params[["methylation_diff"]])
    while (!is.null(dev.list())) dev.off()
  }
      
  ###### Get hypermethylated SMPs/DMRs ######
  for(hh_select in c("hyper","hypo")){
    cat("#####   Calculating ",hh_select,"-methylated ",scope,"s #####\n",sep = "")
    meth_diff[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]] = getMethylDiff(meth_poly[[scope]][[variable]][[context]][[comp_var[2]]],
                                                                                          difference = init_params[["methylation_diff"]],
                                                                                          qvalue = init_params[["alpha"]],
                                                                                          type = hh_select,
                                                                                          save.db = FALSE)
    
  ###### Output Annotations for SMPs/DMRs ######
    cat("#####   Getting annotations from ",hh_select,"-methylated ",scope,"s   #####\n",sep = "")
    if(!isEmpty(meth_diff[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]])){
      # meth_genes[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]] = lapply(init_params[["features"]],function(feature){
      #   res = lapply(X = unique(as.character(meth_diff[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]][["chr"]])),FUN = function(x){
      #   temp = chr_features[[x]][[feature]]
      #   temp2 = temp[temp[["start"]] %in% meth_diff[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]][["start"]],"name"]
      #   return(temp2)
      # })
      #   names(res) = unique(as.character(meth_diff[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]][["chr"]]))
      #   return(res)
      #   })
      # names(meth_genes[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]]) = init_params[["features"]]
      # meth_genes[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]] = meth_genes[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]][!isEmpty(meth_genes[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]])]
      annot = annotateWithGeneParts(as(meth_diff[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]],"GRanges"),
                                    feature = genome_annotation)
      annot_TSS = getAssociationWithTSS(annot)
      meth_genes[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]] = apply(getMembers(annot),MARGIN = 2,function(x){
        unique(annot@dist.to.TSS$feature.name[as.logical(x)])
      })
      names(meth_genes[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]]) = c("promoters","exons","introns")
      
      meth_distTSS[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]] = sapply(colnames(getMembers(annot)),function(x){
        tmp = annot_TSS[as.logical(getMembers(annot)[,x]),,drop = FALSE]
        genes = unique(tmp[,"feature.name"])
        mean_dist = sapply(genes,function(gene){
          cat(x," | ",format(round(100*which(genes %in% gene)/length(genes),digits = 2),nsmall = 2),"%\r",sep = "")
          mean(annot_TSS[annot_TSS[["feature.name"]] %in% gene,"dist.to.feature"])
        })
        cat("\n")
        return(mean_dist)
      })
      names(meth_distTSS[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]]) = c("promoters","exons","introns")
      
      sapply(names(meth_genes[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]]),
             function(feature){
               write.table(x = meth_genes[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]][[feature]],
                           file = paste0("Annotation/",scope,"_",variable,"_",context,"_",comp_var[2],"_",hh_select,"_",feature,".txt"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
             })
      png(filename = paste0("Annotation/DistTSS_",scope,"_",variable,"_",context,"_",comp_var[2],"_",hh_select,".png"),width = 1080,height = 1080,units = "px")
      boxplot(meth_distTSS[[scope]][[variable]][[context]][[comp_var[2]]][[hh_select]],xlab = "Feature",ylab = "Distance to TSS")
      while (!is.null(dev.list())) dev.off()
      
      png(filename = paste0("Annotation/Piechart_",scope,"_",variable,"_",context,"_",comp_var[2],"_",hh_select,".png"),width = 1080,height = 1080,units = "px")
      par(cex = 2.5)
      plotTargetAnnotation(annot)
      while (!is.null(dev.list())) dev.off()
    }
    }
  }
  }
  }
  }
  }
    
  # dir.create("percent_methylation",showWarnings = FALSE)
  # for(i in c("raw","tiles")){
  #   for(f in c("mean","sd")){
  #     for(Annotation_idx in names(perc_meth[[i]][[f]])){
  #       write.table(x = data.frame(perc_meth[[i]][[f]]),file = paste0("./percent_methylation/perc_meth_",context,"_",i,"_",f,"_",Annotation_idx,".txt"),quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
  #     }
  #     perc_PCA[[paste0(i,"_",f)]] = prcomp(perc_meth_list[[paste0("perc_meth_",i,"_",f)]])
  #     perc_PCA[[paste0("var_",i,"_",f)]] = perc_PCA[[paste0(i,"_",f)]]$sdev^2/sum(perc_PCA[[paste0(i,"_",f)]]$sdev^2)
  #     png(filename = paste0("./percent_methylation/PCA_",context,"_",i,"_",f,".png"),width = 1440,height = 810,units = "px")
  #     print(ggplot(as.data.frame(perc_PCA[[paste0(i,"_",f)]]$x),aes(x = PC1, y = PC2, group = perc_meth_list$Annotation)) + 
  #       geom_point(size=4,aes(shape = perc_meth_list$regions,color = perc_meth_list$Annotation)) +
  #       labs(title = paste0("PCA for ",if(i == "raw"){"single"}else{i}," methylation (",f,") for ",context),x=paste0("PC1: ",round(perc_PCA[[paste0("var_",i,"_",f)]][1]*100,1),"%"),y=paste0("PC2: ",round(perc_PCA[[paste0("var_",i,"_",f)]][2]*100,1),"%")) +
  #       theme(legend.position="right") +
  #       scale_color_discrete(name = "Annotation") +
  #       scale_shape_manual(values = c(0,1,15:17),name = "Region") +
  #       stat_ellipse(show.legend = FALSE,geom = "polygon", alpha = 0.25, aes(fill = perc_meth_list$Annotation)))
  #     while (!is.null(dev.list())){dev.off()}
  #   }
  # # }
  # }
  }

save.image(paste0(init_params[["name"]],"_final.RData"))
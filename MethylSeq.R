#######     DNA methylation analysis      ##########

packages=c("methylKit","Rsamtools","genomation","rChoiceDialogs","goseq","BiocParallel","parallel","tools","factoextra",
           "R.utils","ggplot2")

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

options(max.print = 2000000000)
options(stringsAsFactors = FALSE)

###### Cluster commands ######
if(!interactive()){
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  must_args = c("rdata","wd","name","bismark","import","context","t")
  if(!all(must_args %in% names(args))){
    help = matrix(data = c("--rdata    ","RData file path",
                    "--wd","Working directory path",
                    "--name","Experiment name (all data will output to a directory by that name in the working directory)",
                    "--bismark","Bismark output BAM files path",
                    "--import","Import BAM files (Y), do not import BAM files (N) or only import BAM files (O); default = N",
                    "--context","Methylation context: CpG, CHG or CHH",
                    "--savedb","Data will be output into flat files (for reduced RAM usage); default = FALSE",
                    "--t","Number of compute threads",
                    "--arg","Additional R arguments (multiple arguments in separate flags)")
                  ,ncol = 2,byrow = TRUE)
    prmatrix(help,quote = FALSE,rowlab = rep("",nrow(help)),collab = rep("",2))
    stop(paste0("Missing command line input --> ",paste(must_args[!must_args %in% names(args)],collapse = " | ")), call. = TRUE)
  }
  
  cat("Loading RData file: ",args[["rdata"]],"\n", sep = "")
  load(args[["rdata"]])
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  wd = args[["wd"]]
  Experiment_name = args[["name"]]
  section = c("Run dataset on regions (Exons/Introns/Promoters/TSSes)","Run all regions (RAM/CPU intensive)")
  if(args[["import"]] == "O"){
    section = "Import dataset"
  }else{
    if(args[["import"]] == "Y"){
      section =	c(section,"Import dataset")
  }else{
    section =	section[!section %in% "Import dataset"]
    }
  }
  
  dir.create(path = paste0(wd,"/",Experiment_name),showWarnings = FALSE)
  setwd(paste0(wd,"/",Experiment_name))
  cat("Working directory: ",getwd(),"\n", sep = "")
  cat("Experiment name: ",Experiment_name,"\n", sep = "")
  
  experimental_design[,"BAM_file"] = paste0("/",args[["bismark"]],basename(experimental_design[,"BAM_file"]))
  save_db = "savedb" %in% names(args)
  regions = NULL
  if("Run dataset on regions (Exons/Introns/Promoters/TSSes)" %in% section){
    regions = c("promoters","exons","introns","TSSes")
  }
  if("Run all regions (RAM/CPU intensive)" %in% section){
    regions = c(regions,"all_regions")
    }
  context = args[["context"]]
  threads = as.numeric(args[["t"]])
  cat("Number of threads: ",threads,"\n", sep = "")
  
  ##### Additional arguments #####
  add_args = unname(args[which(names(args) %in% "arg")])
  for(i in seq_along(add_args)){
    cat("Performing: ",add_args[[i]],"\n", sep = "")
    eval(parse(text = add_args[[i]]))
  }
  
  ##### Register threads #####
  if(.Platform$OS.type == "unix"){
    register(BPPARAM = MulticoreParam(workers = threads))
  }else{
    register(BPPARAM = SerialParam())
  }
  

##### RData output ######
.classes = NULL
for(.obj in ls()){
  suppressWarnings({.classes[.obj] = class(get(.obj))})
}
prmatrix(matrix(data = c(ls(),.classes),nrow = length(ls()),ncol = 2),quote = FALSE,rowlab = rep("",length(ls())),collab = rep("",2))
rm(.classes,.obj)

}else{
  
  threads = detectCores()
  
  ##### Register threads #####
  if(.Platform$OS.type == "unix"){
    register(BPPARAM = MulticoreParam(workers = threads))
  }else{
    register(BPPARAM = SerialParam())
  }
  
  ###### Set working directory ######
  wd=rchoose.dir(caption = "Choose working directory:")
  setwd(wd)
  
  ##### Setup experiment name ######
  if(!exists("Experiment_name")){Experiment_name = as.character(readline(prompt = "Select experiment name: "))
  if(any(list.files(wd)==Experiment_name)){
    if(rselect.list(choices = c("Overwrite","Create new folder"),multiple = FALSE,title = "Folder exists...")=="Create new folder"){
      Experiment_name_check = Experiment_name
      i=2
      while(any(list.files(wd)==Experiment_name_check)){
        cat("Experiment ",Experiment_name_check," already exists. Changing name...\n", sep = "")
        Experiment_name_check = paste0(Experiment_name,"_(",i,")")
        i=i+1
      }
      Experiment_name=Experiment_name_check
      rm(Experiment_name_check)
    }
  }
  dir.create(Experiment_name,showWarnings = FALSE)
  }
  setwd(paste0(wd,"/",Experiment_name))
  
  section = rselect.list(choices = c("Input settings","Import dataset","Run dataset on regions (Exons/Introns/Promoters/TSSes)","Run all regions (RAM/CPU intensive)"),multiple = TRUE,title = "Select run")
  
  ###### Run settings ######
  if("Input settings" %in% section){
    
    ###### Input analysis settings ######
    cat("#####   Input analysis settings   #####\n")
    if(rselect.list(choices = c("Yes","No"),multiple = FALSE,title = "Save to flat files?") == "Yes"){
      save_db = TRUE
    }else{
      save_db = FALSE
    }
    experimental_design = data.frame(Sample = NA_character_,
                                     Factor = NA_character_,
                                     Factor_name = NA_character_,
                                     BAM_file = rchoose.files(caption = "Select sorted and indexed BAM files: "),
                                     Raw_DB_file = NA_character_,
                                     stringsAsFactors = FALSE)
    assembly = readline(prompt = "Organism: ")
    context = rselect.list(choices = c("CpG","CHG","CHH"),multiple = FALSE,title = "Select contexts to analyze")
    cat("#####   Including following contexts in analysis: ",paste0(context,collapse = ", "),"   #####\n", sep = "")
    window_size=as.numeric(readline(prompt = "Window size for differentially methylated regions (default = 100): "))
    alpha = as.numeric(readline(prompt = "pvalue cutoff for methylation difference (default = 0.05): "))
    methylation_diff=as.numeric(readline(prompt = "Methylation difference cutoff percent (default = 25): "))
    if(is.na(window_size)){window_size = 100}
    if(is.na(alpha)){alpha = 0.05}
    if(is.na(methylation_diff)){methylation_diff = 25}
    
    ###### Input BAM files with metadata (these have to be sorted and indexed with samtools) ######
    cat("#####   Input samples and metadata   #####\n")
    for(i in 1:nrow(experimental_design)){
      experimental_design[[i,"Sample"]] = readline(prompt = paste0("Sample name for ",basename(experimental_design[[i,"BAM_file"]]),": "))
      experimental_design[[i,"Factor"]] = as.numeric(readline(prompt = "Control replicates = 0 / Test replicates = 1,2,3,etc.: "))
    }
    experimental_design = experimental_design[order(experimental_design[,"Factor"]),]
    
    # Sample names
    for(i in as.numeric(unique(experimental_design[,"Factor"]))){
      experimental_design[experimental_design[,"Factor"] == i,"Factor_name"] = readline(prompt = paste0("Factor name for ",paste0(experimental_design[experimental_design[,"Factor"]==i,"Sample"],collapse = " | "),": "))
    }
    cat("#####   Input annotation files   #####\n")
    Annotation_names = rchoose.files(caption = "Choose annotation files: ",multi = TRUE)
    names(Annotation_names) = Annotation_names
    for(i in Annotation_names){
      Annotation_names[i] = readline(prompt = paste0("Annotation name for ",i,": "))
    }
  
    regions = NULL
    if(any(section=="Run dataset on regions (Exons/Introns/Promoters/TSSes)")){regions=rselect.list(choices = c("promoters","exons","introns","TSSes"),multiple = TRUE,title = "Select functional regions to analyze",preselect = c("promoters","exons","introns","TSSes"))}
    if(any(section=="Run all regions (RAM/CPU intensive)")){regions = c(regions,"all_regions")}
    
    ###### Import annotation files ######
    cat("#####   Importing Annotation files   #####\n")
    Annotation_list = list()
    for(i in names(Annotation_names)){
      Annotation_list[[Annotation_names[i]]] = readTranscriptFeatures(location = i)
    }
    rm(Annotation_names)
    save.image(paste0(Experiment_name,"_input.RData"))
  }
}


#####   Run files in selected contexts   ######
  start_time=Sys.time()

#####  Import samples  #####
if("Import dataset" %in% section){
  cat("#####   Importing sample files in ",context," context   #####\n", sep = "")
  if(save_db){
  tryCatch({processBismarkAln(location = as.list(experimental_design[,"BAM_file"]),sample.id = as.list(experimental_design[,"Sample"]),treatment = as.numeric(experimental_design[,"Factor"]),assembly = assembly,save.context = context,read.context = context,save.folder = paste0("./MethylDB/raw_meth_",context,"/"),save.db = TRUE)},error=function(e){cat("ERROR : ",conditionMessage(e), "\n", sep = "")})
  for(i in 1:nrow(experimental_design)){
    experimental_design[[i,"Raw_DB_file"]]=paste0("./MethylDB/raw_meth_",context,"/",experimental_design[[i,"Sample"]],"_",context,".txt")
  }
  write.table(x = data.frame(experimental_design[,c("Sample","Factor")],Factor_name = basename(experimental_design[,"Factor_name"])),file = "./Experimental_design.txt",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
  raw_meth = methRead(location = as.list(experimental_design[,"Raw_DB_file"]),sample.id = as.list(experimental_design[,"Sample"]),treatment = as.numeric(experimental_design[,"Factor"]),assembly = assembly,dbtype = "tabix",dbdir = paste0("./MethylDB/raw_meth_",context,"/"),context = context,pipeline = "bismark")
  gc()
  }else{
  raw_meth = processBismarkAln(location = as.list(experimental_design[,"BAM_file"]),sample.id = as.list(experimental_design[,"Sample"]),treatment = as.numeric(experimental_design[,"Factor"]),assembly = assembly,read.context = context)
  names(raw_meth) = experimental_design[,"Sample"]
  for(i in 1:nrow(experimental_design)){
    experimental_design[[i,"Raw_DB_file"]] = paste0("./MethylDB/raw_meth_",context,"/",experimental_design[[i,"Sample"]],"_",context,".txt")
  }
  
  save.image(paste0(Experiment_name,"_",context,"_imported_data.RData"))
  }
}
  
  if("Run dataset on regions (Exons/Introns/Promoters/TSSes)" %in% section | "Run all regions (RAM/CPU intensive" %in% section){
  ######  Basic stats in numbers  ######
  cat("#####   Calculating basic stats, plots and coverage plots   #####\n")
  dir.create("Methyl_stats",showWarnings = FALSE)
  dir.create("Methyl_stats_plots",showWarnings = FALSE)
  dir.create("Coverage_stats_plots",showWarnings = FALSE)
  for(i in experimental_design[,"Sample"]) {
    sink(paste0("./Methyl_stats/",context,"_getMethylStats.txt"),append = TRUE)
    print(i)
    getMethylationStats(raw_meth[[i]],plot = FALSE,both.strands = FALSE)
    sink()
    
    ###### Basic stats plots ######
    png(filename = paste0("./Methyl_stats_plots/",context,"_Methyl_stats_",i,".png"),width = 1920,height = 1080,units = "px")
    getMethylationStats(raw_meth[[i]],plot = TRUE,both.strands = FALSE)
    while (!is.null(dev.list())){dev.off()}
    
    ###### Coverage stats ######
    png(filename = paste0("./Coverage_stats_plots/",context,"_Coverage_stats_",i,".png"),width = 1920,height = 1080,units = "px")
    getCoverageStats(raw_meth[[i]],plot = TRUE,both.strands = FALSE)
    while (!is.null(dev.list())){dev.off()}
  }
  
  ###### Filter and normalize each sample by read coverage ######
  cat("#####   Filtering data   #####\n")
  raw_meth_filter = filterByCoverage(methylObj = raw_meth,lo.count = 10,lo.perc = NULL,hi.count = NULL,hi.perc = 99.9,save.db = save_db,suffix = paste0(context,"_filtered"),dbdir = paste0("./MethylDB/raw_meth_",context,"/"))
  cat("#####   Normalizing data   #####\n")
  raw_meth_norm = normalizeCoverage(obj = raw_meth_filter, method = "median",save.db = save_db,suffix = paste0(context,"_norm"),dbdir = paste0("./MethylDB/raw_meth_",context,"/"))
  
  ###### Create tiles ######
  cat("#####   Creating region tiles   #####\n")
  tiles = tileMethylCounts(raw_meth_norm,win.size = window_size,step.size = window_size,mc.cores = threads,save.db = save_db,dbdir = paste0("./MethylDB/raw_meth_",context,"/"),suffix = paste0(context,"_tiles"))
  
  ###### Unite samples ######
  cat("#####   Uniting samples   #####\n")
  if("Run dataset on regions (Exons/Introns/Promoters/TSSes)" %in% section){
    if(context=="CpG"){
      raw_meth_unite = unite(object = raw_meth_norm,destrand = TRUE,save.db = FALSE,mc.cores = threads)
      tiles_unite = unite(object = tiles,destrand = TRUE,save.db = FALSE,mc.cores = threads)
    } else {
      raw_meth_unite = unite(object = raw_meth_norm,destrand = FALSE,save.db = FALSE,mc.cores = threads)
      tiles_unite = unite(object = tiles,destrand = FALSE,save.db = FALSE,mc.cores = threads)
    }
  }
  
  ###### Data lists ######
  raw_meth_unite_region = list()
  tiles_unite_region = list()
  perc_meth = list()
  perc_meth_list = list(perc_meth_raw_mean = data.frame(),
                        perc_meth_raw_sd = data.frame(),
                        perc_meth_tiles_mean = data.frame(),
                        perc_meth_tiles_sd = data.frame(),
                        Annotation = rep(names(Annotation_list),each = length(regions)),
                        regions = rep(regions,length(Annotation_list)))
  perc_PCA = list()
  PCA_raw = list()
  PCA_tiles = list()
  subset_list = list()
  SMP = list()
  DMR = list()
  SMP_diff = list()
  DMR_diff = list()
  SMP_genes = list()
  DMR_genes = list()
  
  ###### Annotation and Regions loop ######
  for(Annotation_idx in names(Annotation_list)){
  for(regions_idx in regions){
  cat("######  Running analysis for ",regions_idx," in ",Annotation_idx,"  ######\n", sep = "")
  
    ###### Extract regions ######
    cat("#####   Extracting regions for ",regions_idx," in ",Annotation_idx,"   #####\n", sep = "")
    if(regions_idx != "all_regions"){
      raw_meth_unite_region[[Annotation_idx]][[regions_idx]] = regionCounts(raw_meth_unite,Annotation_list[[Annotation_idx]][[regions_idx]],save.db = save_db,suffix = paste0(context,"_",Annotation_idx,"_",regions_idx),dbdir = "MethylDB/Unite/",mc.cores = threads)
      tiles_unite_region[[Annotation_idx]][[regions_idx]] = regionCounts(tiles_unite,Annotation_list[[Annotation_idx]][[regions_idx]],save.db = save_db,suffix = paste0(context,"_",Annotation_idx,"_",regions_idx,"_tiles"),dbdir = "MethylDB/Unite/",mc.cores = threads)
    }else{
      raw_meth_unite_region[[Annotation_idx]][[regions_idx]] = regionCounts(raw_meth_unite,Annotation_list[[Annotation_idx]]@unlistData,save.db = save_db,suffix = paste0(context,"_",Annotation_idx,"_",regions_idx),dbdir = "MethylDB/Unite/",mc.cores = threads)
      tiles_unite_region[[Annotation_idx]][[regions_idx]] = regionCounts(tiles_unite,Annotation_list[[Annotation_idx]]@unlistData,save.db = save_db,suffix = paste0(context,"_",Annotation_idx,"_",regions_idx,"_tiles"),dbdir = "MethylDB/Unite/",mc.cores = threads)
    }
    
    dir.create("Sequence_regions",showWarnings = FALSE)
    if(length(raw_meth_unite_region[[Annotation_idx]][[regions_idx]]@row.names) > 0){
      temp = data.frame(raw_meth_unite_region[[Annotation_idx]][[regions_idx]]@.Data)
      names(temp) = raw_meth_unite_region[[Annotation_idx]][[regions_idx]]@names
      write.table(x = temp,file = paste0("./Sequence_regions/Single_",Annotation_idx,"_",regions_idx,"_",context,".txt"),quote = FALSE,sep = "\t",row.names = FALSE)
      rm(temp)
      temp = percMethylation(raw_meth_unite_region[[Annotation_idx]][[regions_idx]])
      perc_meth[["raw"]][["mean"]][[Annotation_idx]][[regions_idx]] = colMeans(percMethylation(methylBase.obj = raw_meth_unite_region[[Annotation_idx]][[regions_idx]]))
      perc_meth[["raw"]][["sd"]][[Annotation_idx]][[regions_idx]] = apply(X = percMethylation(methylBase.obj = raw_meth_unite_region[[Annotation_idx]][[regions_idx]]),MARGIN = 2,FUN = sd)
    }
    if(length(tiles_unite_region[[Annotation_idx]][[regions_idx]]@row.names) > 0){  
      temp = data.frame(tiles_unite_region[[Annotation_idx]][[regions_idx]]@.Data)
      names(temp) = tiles_unite_region[[Annotation_idx]][[regions_idx]]@names
      write.table(x = temp,file = paste0("./Sequence_regions/Tiles_",Annotation_idx,"_",regions_idx,"_",context,".txt"),quote = FALSE,sep = "\t",row.names = FALSE)
      rm(temp)
      perc_meth[["tiles"]][["mean"]][[Annotation_idx]][[regions_idx]] = colMeans(percMethylation(methylBase.obj = tiles_unite_region[[Annotation_idx]][[regions_idx]]))
      perc_meth[["tiles"]][["sd"]][[Annotation_idx]][[regions_idx]] = apply(X = percMethylation(methylBase.obj = tiles_unite_region[[Annotation_idx]][[regions_idx]]),MARGIN = 2,FUN = sd)
    }
  
  numCs = NULL
    for(i in grep(pattern = "numCs",raw_meth_unite_region[[Annotation_idx]][[regions_idx]]@names)){
    numCs = c(numCs,sum(raw_meth_unite_region[[Annotation_idx]][[regions_idx]][[i]]))
  }
  if(!any(numCs == 0)){
    
  ###### Correlation ######
  cat("#####   Creating correlation stats and plots   #####\n")
  dir.create("Correlation",showWarnings = FALSE)
  sink(paste0("./Correlation/",context,"_getCorrelation_",Annotation_idx,"_",regions_idx,".txt"))
  getCorrelation(raw_meth_unite_region[[Annotation_idx]][[regions_idx]],plot = FALSE)
  sink()
  png(paste0("./Correlation/",context,"_Correlation_",Annotation_idx,"_",regions_idx,".png"),width = 1920,height = 1080,units = "px")
  getCorrelation(raw_meth_unite_region[[Annotation_idx]][[regions_idx]],plot = TRUE)
  while (!is.null(dev.list())){dev.off()}
  
  sink(paste0("./Correlation/tiles_",context,"_getCorrelation_",Annotation_idx,"_",regions_idx,".txt"))
  getCorrelation(tiles_unite_region[[Annotation_idx]][[regions_idx]],plot = FALSE)
  sink()
  png(filename = paste0("./Correlation/tiles_",context,"_Correlation_",Annotation_idx,"_",regions_idx,".png"),width = 1920,height = 1080,units = "px")
  getCorrelation(tiles_unite_region[[Annotation_idx]][[regions_idx]],plot = TRUE)
  while (!is.null(dev.list())){dev.off()}
  
  
  ###### Cluster the samples ######
  cat("#####   Creating cluster plots   #####\n")
  dir.create("Cluster_Samples",showWarnings = FALSE)
  png(filename = paste0("./Cluster_Samples/",context,"_Cluster_",Annotation_idx,"_",regions_idx,".png"),width = 1080,height = 1080,units = "px")
  clusterSamples(raw_meth_unite_region[[Annotation_idx]][[regions_idx]],dist = "correlation",method = "ward",plot = TRUE)
  while (!is.null(dev.list())){dev.off()}
  
  png(filename = paste0("./Cluster_Samples/tiles_",context,"_Cluster_",Annotation_idx,"_",regions_idx,".png"),width = 1080,height = 1080,units = "px")
  clusterSamples(tiles_unite_region[[Annotation_idx]][[regions_idx]],dist = "correlation",method = "ward",plot = TRUE)
  while (!is.null(dev.list())){dev.off()}
  
  ###### PCA variances and plots ######
  cat("#####   Creating PCA variances and plots   #####\n")
  dir.create("PCA",showWarnings = FALSE)
  png(filename = paste0("./PCA/",context,"_PCA_Variances_",Annotation_idx,"_",regions_idx,".png"),width = 1920,height = 1080,units = "px")
  PCASamples(raw_meth_unite_region[[Annotation_idx]][[regions_idx]], screeplot=TRUE)
  while (!is.null(dev.list())){dev.off()}
  PCA_raw[[Annotation_idx]][[regions_idx]] = PCASamples(raw_meth_unite_region[[Annotation_idx]][[regions_idx]],obj.return = TRUE)
  png(filename = paste0("./PCA/",context,"_PCA_plot_",Annotation_idx,"_",regions_idx,".png"),width = 1920,height = 1080,units = "px")
  print(fviz_pca_ind(X = PCA_raw[[Annotation_idx]][[regions_idx]],repel = TRUE,habillage = experimental_design[,"Factor_name"],title = paste0("PCA for ",context," ",Annotation_idx," ",regions_idx),labelsize = 8, pointsize = 3)+theme(title = element_text(size = 20),axis.title = element_text(size = 20),axis.text = element_text(size = 20),legend.text = element_text(size = 20)))
  while (!is.null(dev.list())){dev.off()}
  
  png(filename = paste0("./PCA/tiles_",context,"_PCA_Variances_",Annotation_idx,"_",regions_idx,".png"),width = 1920,height = 1080,units = "px")
  PCASamples(tiles_unite_region[[Annotation_idx]][[regions_idx]], screeplot=TRUE)
  while (!is.null(dev.list())){dev.off()}
  PCA_tiles[[Annotation_idx]][[regions_idx]] = PCASamples(tiles_unite_region[[Annotation_idx]][[regions_idx]],obj.return = TRUE)
  png(filename = paste0("./PCA/tiles_",context,"_PCA_plot_",Annotation_idx,"_",regions_idx,".png"),width = 1920,height = 1080,units = "px")
  print(fviz_pca_ind(X = PCA_tiles[[Annotation_idx]][[regions_idx]],repel = TRUE,habillage = experimental_design[,"Factor_name"],title = paste0("PCA for tiles ",context," ",Annotation_idx," ",regions_idx),labelsize = 8, pointsize = 3)+theme(title = element_text(size = 20),axis.title = element_text(size = 20),axis.text = element_text(size = 20),legend.text = element_text(size = 20)))
  while (!is.null(dev.list())){dev.off()}
  
  
  ###### Annotation lists ######
  Annotation_df = list()
  if(regions_idx != "all_regions"){
    cat("#####   Creating Annotation list for ",regions_idx," in ",Annotation_idx,"   #####\n", sep = "")
    if(regions_idx == "promoters"){
    Annotation_lengths = c(0,cumsum(Annotation_list[[Annotation_idx]]$TSSes@seqnames@lengths))
    Chr.start.tss = lapply(X = seq_along(Annotation_list[[Annotation_idx]]$TSSes@seqnames@lengths),FUN = function(x){
      data.frame(Start = Annotation_list[[Annotation_idx]]$TSSes@ranges@start[(Annotation_lengths[x]+1):(Annotation_lengths[x+1])], Gene = Annotation_list[[Annotation_idx]]$TSSes$name[(Annotation_lengths[x]+1):(Annotation_lengths[x+1])])
    })
    names(Chr.start.tss) = as.character(Annotation_list[[Annotation_idx]]$TSSes@seqnames@values)
    
    Annotation_lengths = c(0,cumsum(Annotation_list[[Annotation_idx]]$promoters@seqnames@lengths))
    Chr.start.promoters = lapply(X = seq_along(Annotation_list[[Annotation_idx]]$promoters@seqnames@lengths),FUN = function(x){
      data.frame(Start = unique(Annotation_list[[Annotation_idx]]$promoters@ranges@start[(Annotation_lengths[x]+1):(Annotation_lengths[x+1])]))
    })
    names(Chr.start.promoters) = as.character(Annotation_list[[Annotation_idx]]$promoters@seqnames@values)
    
    for(i in names(Chr.start.promoters)){
      cat("\rGetting genes for promoters in chromosome ",i, sep = "")
      Annotation_df[[Annotation_idx]][[regions_idx]][[i]] = bplapply(X = Chr.start.promoters[[i]]$Start,FUN = function(j,x = i){
        genes_temp = Chr.start.tss[[x]]$Gene[Chr.start.tss[[x]]$Start == j+1000]
        df_temp = data.frame(Start = rep(j,length(genes_temp)), Gene = genes_temp)
        return(df_temp)
      })
      Annotation_df_temp = unlist(Annotation_df[[Annotation_idx]][[regions_idx]][[i]])
      Annotation_df[[Annotation_idx]][[regions_idx]][[i]] = data.frame(Start = as.numeric(Annotation_df_temp[grepl(pattern = "Start",x = names(Annotation_df_temp))]), Gene = Annotation_df_temp[grepl(pattern = "Gene",x = names(Annotation_df_temp))])
      cat("\n")
    }
    rm(Annotation_df_temp,Annotation_lengths,Chr.start.promoters,Chr.start.tss)
  }else{
    Annotation_lengths = c(0,Annotation_list[[Annotation_idx]][[regions_idx]]@seqnames@lengths)
    Annotation_df[[Annotation_idx]][[regions_idx]] = bplapply(X = seq_along(Annotation_list[[Annotation_idx]][[regions_idx]]@seqnames@lengths),FUN = function(x){
      data.frame(Start = as.integer(Annotation_list[[Annotation_idx]][[regions_idx]]@ranges@start[(Annotation_lengths[x]+1):(Annotation_lengths[x+1])]), Gene = Annotation_list[[Annotation_idx]][[regions_idx]]$name[(Annotation_lengths[x]+1):(Annotation_lengths[x+1])])
    })
    names(Annotation_df[[Annotation_idx]][[regions_idx]]) = as.character(Annotation_list[[Annotation_idx]][[regions_idx]]@seqnames@values)
  }
  }
  
  
  for(factor_name_idx in unique(experimental_design$Factor_name)[-1]){
  ###### Calculate Single methylation polymorphisms (SMPs) ######
    cat("#####   Subsetting into samples ",unique(experimental_design$Factor_name)[1]," and ",factor_name_idx,"   #####\n", sep = "")
    subset_SMP = reorganize(methylObj = raw_meth_unite_region[[Annotation_idx]][[regions_idx]],sample.ids = c(experimental_design$Sample[experimental_design$Factor == 0],experimental_design$Sample[experimental_design$Factor_name == factor_name_idx]),treatment = as.numeric(c(experimental_design$Factor[experimental_design$Factor==0],experimental_design$Factor[experimental_design$Factor_name == factor_name_idx])))
    subset_DMR = reorganize(methylObj = tiles_unite_region[[Annotation_idx]][[regions_idx]],sample.ids = c(experimental_design$Sample[experimental_design$Factor == 0],experimental_design$Sample[experimental_design$Factor_name == factor_name_idx]),treatment = as.numeric(c(experimental_design$Factor[experimental_design$Factor==0],experimental_design$Factor[experimental_design$Factor_name == factor_name_idx])))

    cat("#####   Calculating single methylation polymorphisms (SMPs)   #####\n")
    tryCatch({SMP[[factor_name_idx]][[Annotation_idx]][[regions_idx]] = calculateDiffMeth(subset_SMP,mc.cores = threads,save.db = FALSE)},
             error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if(exists("subset_SMP")){rm(subset_SMP)}
    
    ###### Calculate differentially methylated regions (DMRs) ######
    cat("#####   Calculating differentially methylated regions (DMRs)   #####\n")
    tryCatch({DMR[[factor_name_idx]][[Annotation_idx]][[regions_idx]] = calculateDiffMeth(subset_DMR,mc.cores = threads,save.db = FALSE)},
            error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if(exists("subset_DMR")){rm(subset_DMR)}
 
    SMP_names = c(names(SMP[[factor_name_idx]][[Annotation_idx]]),names(SMP[[factor_name_idx]]),names(SMP))
    DMR_names = c(names(DMR[[factor_name_idx]][[Annotation_idx]]),names(DMR[[factor_name_idx]]),names(DMR))
    
    if(all(c(factor_name_idx,Annotation_idx,regions_idx) %in% SMP_names) & all(c(factor_name_idx,Annotation_idx,regions_idx) %in% DMR_names)){
  ###### Get differential methylation stats ######
  cat("#####   Creating stats and plots for SMPs and DMRs for chromosomes   #####\n")
  dir.create("Differential_methylation",showWarnings = FALSE)
  if(sum((SMP[[factor_name_idx]][[Annotation_idx]][[regions_idx]][["qvalue"]] <= alpha)*(abs(SMP[[factor_name_idx]][[Annotation_idx]][[regions_idx]][["meth.diff"]]) >= methylation_diff)) > 0){
    tryCatch({
      sink(paste0("./Differential_methylation/Diff_meth_per_chr_SMP_",context,"_",regions_idx,"_",Annotation_idx,"_factor_",factor_name_idx,".txt"))
      print(diffMethPerChr(SMP[[factor_name_idx]][[Annotation_idx]][[regions_idx]],plot=FALSE,qvalue.cutoff=alpha, meth.cutoff=methylation_diff))
      sink()},
      error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    png(filename = paste0("./Differential_methylation/Diff_meth_per_chr_SMP_",context,"_",regions_idx,"_",Annotation_idx,"_factor_",factor_name_idx,".png"),width = 1920,height = 1080,units = "px")
    tryCatch({
        diffMethPerChr(SMP[[factor_name_idx]][[Annotation_idx]][[regions_idx]],plot=TRUE,qvalue.cutoff=alpha, meth.cutoff=methylation_diff)},
             error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    while (!is.null(dev.list())){dev.off()}
  }
  if(sum((DMR[[factor_name_idx]][[Annotation_idx]][[regions_idx]][["qvalue"]] <= alpha)*(abs(DMR[[factor_name_idx]][[Annotation_idx]][[regions_idx]][["meth.diff"]]) >= methylation_diff)) > 0){
    tryCatch({
      sink(paste0("./Differential_methylation/Diff_meth_per_chr_DMR_",context,"_",regions_idx,"_",Annotation_idx,"_factor_",factor_name_idx,".txt"))
      print(diffMethPerChr(DMR[[factor_name_idx]][[Annotation_idx]][[regions_idx]],plot=FALSE,qvalue.cutoff=alpha, meth.cutoff=methylation_diff))
      sink()},
      error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  png(filename = paste0("./Differential_methylation/Diff_meth_per_chr_DMR_",context,"_",regions_idx,"_",Annotation_idx,"_factor_",factor_name_idx,".png"),width = 1920,height = 1080,units = "px")
  tryCatch({
    diffMethPerChr(DMR[[factor_name_idx]][[Annotation_idx]][[regions_idx]],plot=TRUE,qvalue.cutoff=alpha, meth.cutoff=methylation_diff)},
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  while (!is.null(dev.list())){dev.off()}
  }
  
  ###### Get hypermethylated SMPs/DMRs ######
  cat("#####   Calculating hyper/hypo-methylated SMPs and DMRs   #####\n")
  for(hh_select in c("hyper","hypo")){
    SMP_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]] = getMethylDiff(SMP[[factor_name_idx]][[Annotation_idx]][[regions_idx]],difference = methylation_diff,qvalue = alpha,type = hh_select,save.db = FALSE)
    DMR_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]] = getMethylDiff(DMR[[factor_name_idx]][[Annotation_idx]][[regions_idx]],difference = methylation_diff,qvalue = alpha,type = hh_select,save.db = FALSE)
 
  ###### Output Annotations for SMPs/DMRs ######
  dir.create("Annotation",showWarnings = FALSE)
    cat("#####   Getting annotations from hyper/hypo-methylated SMPs and DMRs   #####\n")
    if(regions_idx != "all_regions"){
    if(!all(isEmpty(SMP_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]]))){
      SMP_genes[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]] = lapply(X = unique(as.character(SMP_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]][["chr"]])),FUN = function(x){
        Annotation_df[[Annotation_idx]][[regions_idx]][[x]]$Gene[Annotation_df[[Annotation_idx]][[regions_idx]][[x]]$Start %in% SMP_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]][["start"]]]
      })
      names(SMP_genes[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]]) = unique(as.character(SMP_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]][["chr"]]))
      SMP_genes[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]] = SMP_genes[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]][!isEmpty(SMP_genes[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]])]
      write.table(x = unique(unlist(SMP_genes[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]],use.names = FALSE)),file = paste0("./Annotation/SMP_",context,"_",factor_name_idx,"_",Annotation_idx,"_",regions_idx,"_",hh_select,".txt"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
    }
    
    if(!all(isEmpty(DMR_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]]))){
      DMR_genes[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]] = lapply(X = unique(as.character(DMR_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]][["chr"]])),FUN = function(x){
        Annotation_df[[Annotation_idx]][[regions_idx]][[x]]$Gene[Annotation_df[[Annotation_idx]][[regions_idx]][[x]]$Start %in% DMR_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]][["start"]]]
      })
      names(DMR_genes[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]]) = unique(as.character(DMR_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]][["chr"]]))
      DMR_genes[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]] = DMR_genes[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]][!isEmpty(DMR_genes[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]])]
      write.table(x = unique(unlist(DMR_genes[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]],use.names = FALSE)),file = paste0("./Annotation/DMR_",context,"_",factor_name_idx,"_",Annotation_idx,"_",regions_idx,"_",hh_select,".txt"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
    }
    }else{
    ##### Annotation #####
    cat("#####   Annotating differentially methylated regions   #####\n")
    if(!all(isEmpty(SMP_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]]))){
      diffAnn_SMP = annotateWithGeneParts(as(SMP_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]],"GRanges"),Annotation_list[[Annotation_idx]])

      sink(paste0("./Annotation/Annotation_stats_SMP_",hh_select,"_",context,"_",Annotation_idx,"_factor_",factor_name_idx,".txt"))
      print(getTargetAnnotationStats(diffAnn_SMP,percentage=TRUE,precedence=TRUE))
      sink()
    
      png(filename = paste0("./Annotation/Annotation_plot_SMP_",hh_select,"_",context,"_",Annotation_idx,"_factor_",factor_name_idx,".png"),width = 1920,height = 1080,units = "px")
      plotTargetAnnotation(diffAnn_SMP,precedence=TRUE,main=paste0("SMP ",context," (",hh_select," methylation in ",Annotation_idx," database)"))
      while (!is.null(dev.list())){dev.off()}
    }
      
    if(!all(isEmpty(DMR_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]]))){
      diffAnn_DMR = annotateWithGeneParts(as(DMR_diff[[factor_name_idx]][[Annotation_idx]][[regions_idx]][[hh_select]],"GRanges"),Annotation_list[[Annotation_idx]])
      
      sink(paste0("./Annotation/Annotation_stats_DMR_",hh_select,"_",context,"_",Annotation_idx,"_factor_",factor_name_idx,".txt"))
      print(getTargetAnnotationStats(diffAnn_DMR,percentage=TRUE,precedence=TRUE))
      sink()
    
      png(filename = paste0("./Annotation/Annotation_plot_DMR_",hh_select,"_",context,"_",Annotation_idx,"_factor_",factor_name_idx,".png"),width = 1920,height = 1080,units = "px")
      plotTargetAnnotation(diffAnn_SMP,precedence=TRUE,main=paste0("DMR ",context," (",hh_select," methylation in ",Annotation_idx," database)"))
      while (!is.null(dev.list())){dev.off()}
    }
  }
  }
  }
  }
  }
  }
    perc_meth_list$perc_meth_raw_mean = rbind(perc_meth_list$perc_meth_raw_mean,data.frame(matrix(unlist(perc_meth$raw$mean[[Annotation_idx]]),nrow = length(regions),dimnames = list(paste0(Annotation_idx,"_",regions),experimental_design$Sample))))
    perc_meth_list$perc_meth_raw_sd = rbind(perc_meth_list$perc_meth_raw_sd,data.frame(matrix(unlist(perc_meth$raw$sd[[Annotation_idx]]),nrow = length(regions),dimnames = list(paste0(Annotation_idx,"_",regions),experimental_design$Sample))))
    perc_meth_list$perc_meth_tiles_mean = rbind(perc_meth_list$perc_meth_tiles_mean,data.frame(matrix(unlist(perc_meth$tiles$mean[[Annotation_idx]]),nrow = length(regions),dimnames = list(paste0(Annotation_idx,"_",regions),experimental_design$Sample))))
    perc_meth_list$perc_meth_tiles_sd = rbind(perc_meth_list$perc_meth_tiles_sd,data.frame(matrix(unlist(perc_meth$tiles$sd[[Annotation_idx]]),nrow = length(regions),dimnames = list(paste0(Annotation_idx,"_",regions),experimental_design$Sample))))
  }
  for(factor_name_idx in unique(experimental_design$Factor_name)[-1]){
    write.table(unique(c(unlist(SMP_genes[grep(pattern = factor_name_idx,x = names(SMP_genes))]),unlist(DMR_genes[grep(pattern = factor_name_idx,x = names(DMR_genes))]))),file = paste0("./Annotation/Total_diff_meth_genes_",context,"_",factor_name_idx,".txt"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")
  }
  dir.create("percent_methylation",showWarnings = FALSE)
  for(i in c("raw","tiles")){
    for(f in c("mean","sd")){
      for(Annotation_idx in names(perc_meth[[i]][[f]])){
        write.table(x = data.frame(perc_meth[[i]][[f]][[Annotation_idx]]),file = paste0("./percent_methylation/perc_meth_",context,"_",i,"_",f,"_",Annotation_idx,".txt"),quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
      }
      perc_PCA[[paste0(i,"_",f)]] = prcomp(perc_meth_list[[paste0("perc_meth_",i,"_",f)]])
      perc_PCA[[paste0("var_",i,"_",f)]] = perc_PCA[[paste0(i,"_",f)]]$sdev^2/sum(perc_PCA[[paste0(i,"_",f)]]$sdev^2)
      png(filename = paste0("./percent_methylation/PCA_",context,"_",i,"_",f,".png"),width = 1440,height = 810,units = "px")
      print(ggplot(as.data.frame(perc_PCA[[paste0(i,"_",f)]]$x),aes(x = PC1, y = PC2, group = perc_meth_list$Annotation)) + 
        geom_point(size=4,aes(shape = perc_meth_list$regions,color = perc_meth_list$Annotation)) +
        labs(title = paste0("PCA for ",if(i == "raw"){"single"}else{i}," methylation (",f,") for ",context),x=paste0("PC1: ",round(perc_PCA[[paste0("var_",i,"_",f)]][1]*100,1),"%"),y=paste0("PC2: ",round(perc_PCA[[paste0("var_",i,"_",f)]][2]*100,1),"%")) +
        theme(legend.position="right") +
        scale_color_discrete(name = "Annotation") +
        scale_shape_manual(values = c(0,1,15:17),name = "Region") +
        stat_ellipse(show.legend = FALSE,geom = "polygon", alpha = 0.25, aes(fill = perc_meth_list$Annotation)))
      while (!is.null(dev.list())){dev.off()}
    }
  }
  }
  
  
if("Import dataset" %in% section | "Run dataset on regions (Exons/Introns/Promoters/TSSes)" %in% section | "Run all regions (RAM/CPU intensive" %in% section){
  ###### Save and finish ######
  save(list = ls(pattern = "^raw|^tiles"),file = paste0(Experiment_name,"_",context,"_raw_data.RData.bz2"),compress = "bzip2")
  save(list = ls()[!ls() %in% ls(pattern = "^raw|^tiles")],file = paste0(Experiment_name,"_",context,"_results_data.RData"))
  end_time=Sys.time()
  write.table(x = paste0(context,"\t",format(round(end_time-start_time,2),nsmall=2)),file = "./Runtime.txt",append = TRUE,quote = FALSE,row.names = FALSE,col.names = FALSE)
  cat("#####   ",context," context run finished in ",format(round(end_time-start_time,2),nsmall=2),"   #####\n", sep = "")
}
  # End main for-loop
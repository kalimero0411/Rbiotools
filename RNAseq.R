#####    RNA-seq analysis (RSEM / counts)    #####

packages = c("DESeq2","ggplot2","ggrepel","gplots","RColorBrewer","BiocParallel","tximport","readr",
             "pheatmap","goseq","rstudioapi","ReportingTools","factoextra","vegan","rgl","ape","cluster","data.table",
             "parallel","doParallel","RCurl","devtools","GenomicFeatures","apeglm","R.utils","VennDiagram","wordcloud",
             "tm","topGO","Rgraphviz","NOISeq")
if(interactive()){
  packages = c(packages,"rChoiceDialogs")
}
# if(!"bcbioRNASeq" %in% rownames(installed.packages())){
#   install.packages(
#     pkgs = "bcbioRNASeq",
#     repos = c("https://r.acidgenomics.com",getOption("repos"))
#   )
# }
packages = c(packages,"bcbioRNASeq")

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

###### Cluster commands ######
if(!interactive()){
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  must_args = c("rdata","wd","name","process","t")
  if(!all(must_args %in% names(args))){
    help = matrix(data = c("--rdata    ","RData file path",
                           "--wd","Working directory path",
                           "--name","Experiment name (all data will output to a directory",
                           "","by that name in the working directory)",
                           "--input","Aligner output files path (RSEM, STAR or kallisto)",
                           "--design","DESeq2 design formula (e.g. ~Genotype+treatment)",
                           "--reduced","DESeq2 reduced design formula (Default = ~1)",
                           "--process","Process: ",
                           "","1 = Run settings (interactive only)",
                           "","2 = DESeq2",
                           "","3 = LRT",
                           "","4 = Transformation (not parallelized)",
                           "","5 = Heatmaps",
                           "","6 = Dispersion estimates, PCAs and PCoAs",
                           "","7 = Optimize K-means",
                           "","8 = K-means clustering",
                           "","9 = Hierarchial clustering",
                           "","10 = MA-plot",
                           "","11 = DEG heatmap",
                           "","12 = goseq GO analysis",
                           "","13 = topGO analysis",
                           "","14 = Word cloud",
                           "","15 = Venn diagram",
                           "","16 = Variable heatmap and report",
                           "","(e.g. 2,5,6,7,8,9,10,11,12)",
                           "--remove_isoforms","Remove isoform suffix from gene IDs",
                           "--k","Number of clusters for K-means clustering",
                           "--seed","Seed value for random number generator",
                           "--heatmap_no_clust","Cluster rows in heatmaps",
                           "--GO_file","Path to GO annotation file",
                           "--NOISeq","Perform NOISeq correction (ARSyNseq: counts | proportion | FDR)",
                           "--t","Number of compute threads",
                           "--arg","Additional R arguments (multiple arguments in separate flags)")
                  ,ncol = 2,byrow = TRUE)
    prmatrix(help,quote = FALSE,rowlab = rep("",nrow(help)),collab = rep("",2))
    stop(paste0("Missing command line input --> ",paste(must_args[!must_args %in% names(args)],collapse = " | ")), call. = TRUE)
    if(sum(as.numeric(unlist(strsplit(args[["process"]],split = ","))) %in% c(2,3,4)) > 1){
      stop(paste0("Select only one DESeq2 analysis (2,3 or 4)"),call. = TRUE)
    }
  }
  cat("Loading RData file: ",args[["rdata"]],"\n",sep = "")
  load(args[["rdata"]])
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  wd = args[["wd"]]
  Experiment_name = args[["name"]]
  if("input" %in% names(args)){
    input_path = args[["input"]]
  }
  heatmap_row_clust = !"heatmap_no_clust" %in% names(args)
  NOISeq_correction = "NOISeq" %in% names(args)
  dir.create(path = paste0(wd,"/",Experiment_name),showWarnings = FALSE)
  setwd(paste0(wd,"/",Experiment_name))
  if(!exists("design_formula")){
    design_formula = formula(args[["design"]])
  }
  if(!exists("reduced_formula")){
    if("reduced" %in% names(args)){
      reduced_formula = formula(args[["reduced"]])
    }else{
      reduced_formula = formula(~1)
    }
  }
  section = c("Run settings",
              "DESeq2",
              "LRT",
              "Transformation",
              "Heatmaps",
              "Dispersion estimates, PCAs and PCoAs",
              "Optimize K-means",
              "K-means clustering",
              "Hierarchial clustering",
              "MA-plot",
              "DEG heatmap",
              "goseq GO analysis",
              "topGO analysis",
              "Wordcloud",
              "Venn diagram",
              "Variable heatmap and report")
  section = section[as.numeric(unlist(strsplit(args[["process"]],split = ",")))]
  if(!"Optimize K-means" %in% section & "K-means clustering" %in% section){
    if(!"k" %in% names(args)){
      stop(paste0("Error: Please specify the number of k clusters with --k OR K cluster optimization with --process 8"),call. = TRUE)
    }
    k_clusters = as.numeric(args[["k"]])
    cat("Number of K clusters: ",k_clusters,"\n",sep = "")
  }else{
    cat("Performing K cluster optimization\n",sep = "")
  }
  if("GO_file" %in% names(args)){
    GO_file = args[["GO_file"]]
  }
  if("remove_isoforms" %in% names(args)){
    remove_isoforms = TRUE
    isoform_pattern = unname(args[["remove_isoforms"]])
  }else{
    remove_isoforms = FALSE
  }
  if("seed" %in% names(args)){
    random_seed = as.numeric(args[["seed"]])
  }else{
    random_seed = 123
  }
  threads = as.numeric(args[["t"]])
  
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
  cat("Running sections: ",paste(section,collapse = " | "),"\n",sep = "")
  
  ##### Register threads #####
  if(.Platform$OS.type == "unix"){
    register(BPPARAM = MulticoreParam(workers = threads))
  }else{
    register(BPPARAM = SnowParam(workers = threads))
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
  register(BPPARAM = SnowParam(workers = threads))
}

###### Set working directory ######
cat("#####   Select working directory   #####\n")
wd = rstudioapi::selectDirectory(caption = "Choose working directory:")
setwd(wd)

##### Setup experiment name ######
if(!exists("Experiment_name")){Experiment_name = as.character(readline(prompt = "Select experiment name: "))
if(any(list.files(wd) == Experiment_name)){
  if(select.list(choices = c("Overwrite","Create new folder"),multiple = FALSE,title = "Folder exists...",graphics = TRUE) == "Create new folder"){
    Experiment_name_check = Experiment_name
    i = 2
    while(any(list.files(wd) == Experiment_name_check)){
      cat("Experiment ",Experiment_name_check," already exists. Changing name...\n",sep = "")
      Experiment_name_check = paste0(Experiment_name,"_(",i,")")
      i = i+1
    }
    Experiment_name = Experiment_name_check
    rm(Experiment_name_check)
  }
}
}
dir.create(Experiment_name,showWarnings = FALSE)
setwd(paste0(wd,"/",Experiment_name))

section = select.list(choices = c("Run settings",
                                  "DESeq2",
                                  "LRT",
                                  "Transformation",
                                  "Heatmaps",
                                  "Dispersion estimates, PCAs and PCoAs",
                                  "Optimize K-means",
                                  "K-means clustering",
                                  "Hierarchial clustering",
                                  "MA-plot",
                                  "DEG heatmap",
                                  "goseq GO analysis",
                                  "topGO analysis",
                                  "Wordcloud",
                                  "Venn diagram",
                                  "Variable heatmap and report"),
                       multiple = TRUE,
                       title = "Selection sections",
                       graphics = TRUE)

#####    Run settings    #####
if("Run settings" %in% section){
  # Get mapper output and factor number
  Mapper = select.list(choices = c("RSEM","Kallisto","HTseq-count","Counts"),multiple = FALSE,title = "Select Mapper",graphics = TRUE)
  input_path = rstudioapi::selectDirectory(caption = "Choose mapping output directory: ")
  
  # Lengths from GTF file
  if(Mapper == "Counts" | Mapper == "HTseq-count"){
    genes_isoforms = "genes"
    if(select.list(choices = c("Yes","No"),multiple = FALSE,title = "Lengths from gtf file?",graphics = TRUE) == "Yes"){
    gtf = read.table(file = file.choose(),header = FALSE,sep = "\t",comment.char = "#")
    gtf = gtf[grep(pattern = "transcript|gene",gtf[,3],ignore.case = TRUE),]
    if(select.list(choices = c("Genes","Transcripts"),multiple = FALSE,title = "Genes or transcript IDs?",graphics = TRUE) == "Genes"){
      genes_transcripts = "gene_id"
    }else{
      genes_transcripts = "transcript_id"
    }
    gtf[["names"]] = sapply(1:nrow(gtf), function(x){
        vec = sub(pattern = ";",replacement = "",strsplit(gtf[x,9],split = "; ")[[1]])
        vec = sub(vec[grep(pattern = genes_transcripts,vec)],pattern = paste0(genes_transcripts," "),replacement = "")
        })
    gtf = gtf[!duplicated(gtf[["names"]]),]
    gene_length = data.frame(length = gtf[,5] - gtf[,4],row.names = gtf[["names"]])
  }
    }
  if(select.list(choices = c("Yes","No"),multiple = FALSE,title = "Remove isoforms for annotation?",graphics = TRUE) == "Yes"){
    remove_isoforms = TRUE
    isoform_pattern = readline(prompt = "Input regular expression for isoform removal: ")
  }else{
    remove_isoforms = FALSE
  }
  
  if(Mapper == "Kallisto"){
    genes_isoforms = "genes"
    tx2gene = read.table(selectFile(caption = "Select transcript to gene file",path = getwd()))
  }
  
  if(Mapper == "Counts"){
    count_genes = as.numeric(readline(prompt = "Select gene column: "))
    if(!exists("gtf")){
      count_length = as.numeric(readline(prompt = "Select length column: "))
    }
    count_column = as.numeric(eval(expr = parse(text = readline(prompt = "Select counts column (set range e.g. 4:7): "))))
  }
  
  # Create experimental design
  cat("Text file: column 1 sample names with title, column 2 file path with title, and the first row of the rest of the columns as factor names\n")
  meta_input = select.list(choices = c("Text file","Manual input"),multiple = FALSE,title = "Input method",graphics = TRUE)
  if(meta_input == "Text file"){
    experimental_design = read.table(file = file.choose(),header = TRUE,sep = "\t",row.names = 1)
    if(!all(experimental_design[,1] %in% list.files(path = input_path,full.names = TRUE))){
      stop("Error: file names in experimental design file must match file names (multiple files) / sample names (one file counts)\n",call. = TRUE)
    }
    file_names = experimental_design[,1]
    experimental_design = experimental_design[,-1,drop = FALSE]
    factors = colnames(experimental_design)
    for(i in colnames(experimental_design)){
      experimental_design[[i]] = as.factor(experimental_design[[i]])
    }
    if(Mapper == "RSEM"){
      genes_isoforms = if(grepl(pattern = "[.]genes[.]results",file_names[1])){
        "genes"
      }else{
        "isoforms"
      }
    }
  }else{
    
    # RSEM input
    if(Mapper == "RSEM"){
      genes_isoforms = select.list(choices = c("genes","isoforms"),multiple = FALSE,title = "Select analysis of genes or isoforms",graphics = TRUE)
      file_names = list.files(path = input_path,pattern = paste0(".",genes_isoforms,".results"),full.names = TRUE)
    }
    
    # Kallisto input
    if(Mapper == "Kallisto"){
      file_names = list.files(path = input_path,pattern = "_abundance.tsv",full.names = TRUE)
    }
    
    if(Mapper == "Counts" | Mapper == "HTseq-count"){
      file_names = rchoose.files(caption = "Select count file(s): ",multi = TRUE)
      
      # HTseq-count input
      # if(Mapper == "HTseq-count"){
      #   if(Sys.info()[["sysname"]] == "Linux"){
      #     for(i in file_names){
      #       try(system(command = paste0("sed -i -E 's/^N_.{1,}$//g' ",i," && sed -i -E '/^$/d' ",i)))
      #     }
      #   }else{
      #     cat("#################################################################\n")
      #     cat("Make sure to remove the first 4 lines of the output if they exist\n")
      #     cat("N_unmapped | N_multimapping | N_noFeature | N_ambiguous\n")
      #     cat("#################################################################\n")
      #     invisible(readline(prompt="Press [enter] to continue"))
      #   }
      # }
      
    }
    
  factors = c()
  while(TRUE){
  factors = c(factors,readline(prompt = "Add factor [blank = finish]: "))
  if(factors[length(factors)] == ""){
    factors = factors[-length(factors)]
    break
  }
  }
  experimental_design = as.data.frame(matrix(data = NA,
                                             nrow = length(file_names),
                                             ncol = length(factors),
                                             dimnames = list(file_names,factors)))
  
  # Select sample names
  if(select.list(choices = c("Yes","No"),multiple = FALSE,title = "Select sample names?",graphics = TRUE)=="Yes"){
    for(i in 1:nrow(experimental_design)){
      rownames(experimental_design)[i] = readline(prompt = paste0("Select name for sample ",rownames(experimental_design)[i],": "))
    }
  }
  
  # Insert factors
  for(i in factors){
    for(j in 1:nrow(experimental_design)){
      experimental_design[j,i] = readline(prompt = paste0("Input ",i," for sample ",rownames(experimental_design)[j],": "))
    }
    experimental_design[[i]] = as.factor(experimental_design[[i]])
  }
  }
  
  # Select variable colors
  color_select = list()
  if(select.list(choices = c("Yes","No"),multiple = FALSE,title = "Select colors for variables?",graphics = TRUE)=="Yes"){
    color_list = NULL
    for(f in colnames(experimental_design)){
      for(i in unique(experimental_design[,f])){
        color_list[i] = select.list(choices = colors(),multiple = FALSE,title = paste0("Select color for variable ",i,": "),graphics = TRUE)
      }
      color_select[[f]] = color_list
    }
    rm(color_list)
  }
  
  # Select variable reference
  primary_factor = select.list(choices = colnames(experimental_design),multiple = FALSE,title = "Select primary factor",graphics = TRUE)
  variable_ref = c()
  for(f in colnames(experimental_design)){
    variable_ref[f] = select.list(choices = as.character(unique(experimental_design[,f])),multiple = FALSE,title = paste0("Select control variable"),graphics = TRUE)
  }
  experimental_design = cbind(experimental_design[,colnames(experimental_design)!=primary_factor,drop = FALSE],experimental_design[,primary_factor,drop = FALSE])
  write.table(x = cbind(Sample_name = rownames(experimental_design),File_name = file_names,experimental_design),file = "Experimental_design.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  cat("Independent design: ~",paste(factors,collapse = " + "),"\n",sep = "")
  cat("Interaction term = Factor1:Factor2 \n",sep = "")
  design_formula = formula(readline(prompt = "Select main design formula: "))
  cat("Reduced design, such as batch factor. No design = ~1","\n",sep = "")
  reduced_formula = formula(readline(prompt = "Select reduced design formula: "))
  if(select.list(choices = c("Yes","No"),multiple = FALSE,title = "Cluster rows in heatmaps?",graphics = TRUE)=="Yes"){
    heatmap_row_clust = TRUE
  }else{
    heatmap_row_clust = FALSE
    }
alpha = as.numeric(readline(prompt = "Choose alpha cutoff: "))
FDR_cutoff = as.numeric(readline(prompt = "Choose FDR cutoff: "))
rlog_vst = select.list(choices = c("rlog","VST"),multiple = FALSE,title = "Select transformation",graphics = TRUE)
NOISeq_correction = select.list(choices = c("Yes","No"),multiple = FALSE,title = "NOISeq correction?",graphics = TRUE) == "Yes"
random_seed = as.numeric(readline(prompt = "Select seed value for random number generator: "))

if(all(!"Optimize K-means" %in% section,"K-means clustering" %in% section)){
  k_clusters = as.numeric(readline(prompt = "Select number of clusters (k): "))
  }

if(select.list(choices = c("Yes","No"),multiple = FALSE,title = "Load GO annotations?",graphics = TRUE) == "Yes"){
  GO_file = rchoose.files(caption = "Choose GO annotation file:")
}
    save.image(paste0(Experiment_name,".RData"))
}
  }

if(any(!section %in% "Run settings")){
##### Load GO annotation #####
if(exists("GO_file")){
  cat("##### Downloading GO terms #####\n")
  GO_terms = c(as.list(GO.db::GOTERM),as.list(GO.db::GOSYNONYM))
  delim_fun = function(GO_file){
    GO_temp = read.table(GO_file,sep = "\t")
    if(any(lengths(gregexpr(pattern = "GO",GO_temp[,2])) > 1)){
      GO_first = GO_temp[lengths(gregexpr(pattern = "GO",GO_temp[,2])) > 1,2][1]
      return(readMappings(GO_file,IDsep = substr(x = GO_first,start = 11,stop = unlist(gregexpr(pattern = "GO:[0-9]{7}",GO_first))[2]-1)))
    }else{
      geneGO_temp = readMappings(GO_file)
      return(split(unlist(geneGO_temp, use.names = FALSE), names(geneGO_temp)))
    }
  }
  geneGO = delim_fun(GO_file)
  GO_table = data.frame(Gene = rep(names(lengths(geneGO)),lengths(geneGO)),GO_ID = unlist(geneGO,use.names = FALSE))
  # GO_table = data.frame(Gene = unique(GO_table$V1), GO_ID = unlist(bplapply(unique(GO_table$V1),function(x){
  #   return(paste(GO_table[GO_table$V1 %in% x,"V2"],collapse = ";"))
  # })))
  cat("##### Applying GO terms #####\n")
  GO_table_temp = as.data.frame(matrix(unlist(bplapply(unique(GO_table[["GO_ID"]]),function(x){
    if(x %in% names(GO_terms)){
      return(c(attr(GO_terms[[x]],which = "Term",exact = TRUE),
               attr(GO_terms[[x]],which = "Ontology",exact = TRUE),
               attr(GO_terms[[x]],which = "Definition",exact = TRUE)))
    }else{
      return(rep(NA_character_,3))
    }})),
    ncol = 3,byrow = TRUE),row.names = unique(GO_table[["GO_ID"]]))
  GO_table = cbind(GO_table,GO_table_temp[match(GO_table[["GO_ID"]],rownames(GO_table_temp)),])
  rm(GO_table_temp)
  colnames(GO_table)[3:5] = c("Term","Ontology","Definition")
  write.table(x = GO_table,file = "GO_table_info.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
  write.table(x = GO_table[,c("Gene","GO_ID")],file = "GO_table.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
}

##### Load custom functions #####
eval(parse(text = getURL("https://raw.githubusercontent.com/kevinblighe/clusGapKB/master/clusGapKB.R", ssl.verifypeer = FALSE)))

# PCA plot function for PC2 and PC3
plotPCA_PC123 = function (object, intgroup = "condition", ntop = 500,returnData = FALSE){
  rv = rowVars(assay(object))
  if(length(rv) < 3){
    cat("#### Not enough genes to create PCA for: ",intgroup," ####\n")
    return(NULL)
  }
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  pca = prcomp(t(assay(object)[select, ]))
  percentVar = pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))){
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df = as.data.frame(colData(object)[, intgroup,drop = FALSE])
  group = if (length(intgroup) > 1){
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }else{
    colData(object)[[intgroup]]
  }
  d = data.frame(PC1 = pca$x[, 1],PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group,intgroup.df, name = colnames(object))
  if (returnData){
    attr(d, "percentVar") = percentVar[1:3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group")) + geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[2]*100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[3]*100), "% variance")) + coord_fixed()
}

# 3D PCA plot function
PCA_3D = function(PC_factor,pca_type){
  multi_done = list()
  for(run_factor in PC_factor){
    colors_3d = brewer.pal(n = 9,name = "Set1")[match(PCA_data[["Single_factor"]][[run_factor]][["group"]],table = unique(PCA_data[["Single_factor"]][[run_factor]][["group"]]))]
    # Single factor 3D PCAs
    if("Single_factor" %in% names(PCA_data) & grepl(pattern = 1,pca_type)){
      dir.create("animation_merge",showWarnings = FALSE)
      for(degree in 1:360) {
        open3d()
        par3d(windowRect = c(20, 30, 1080, 1080),dev = unname(rgl.dev.list()))
        plot3d(x = PCA_data[["Single_factor"]][[run_factor]][["PC1"]],
               y = PCA_data[["Single_factor"]][[run_factor]][["PC2"]],
               z = PCA_data[["Single_factor"]][[run_factor]][["PC3"]],
               col = colors_3d,
               size = 2,
               xlab = paste0("PC1: ",round(100 * attr(PCA_data[["Single_factor"]][[run_factor]], "percentVar"))[1],"% variance"),
               ylab = paste0("PC2: ",round(100 * attr(PCA_data[["Single_factor"]][[run_factor]], "percentVar"))[2],"% variance"),
               zlab = paste0("PC3: ",round(100 * attr(PCA_data[["Single_factor"]][[run_factor]], "percentVar"))[3],"% variance"),
               type = "s") +
          legend3d("topright", legend = unique(PCA_data[["Single_factor"]][[run_factor]][["group"]]), col = brewer.pal(n = 9,name = "Set1")[1:length(unique(PCA_data[["Single_factor"]][[run_factor]][["group"]]))], pch = 16, cex=1, inset=c(0.02))
        view3d(userMatrix=rotationMatrix(2*pi * degree/360, 0, 1, 0))
        rgl.snapshot(filename=paste("animation_merge/frame-",
                                    sprintf("%03d", degree), ".png", sep=""))
        while(length(rgl.dev.list()) != 0){rgl.close()}
        cat("\rRunning single factor ",run_factor," | ",format(round((degree/360)*100,digits = 2),nsmall = 2),"%", sep = "")
      }
      cat("\nRunning ffmpeg...\n")
      try(system(paste0("ffmpeg -y -hide_banner -loglevel warning -r 60 -y -i ",getwd(),"/animation_merge/frame-%03d.png ./",rlog_vst,"/PCA/PCA_3D_",rlog_vst,"_",genes_isoforms,"_",run_factor,".mp4")))
      unlink("animation_merge",recursive = TRUE)
    }
    
    # Multiple factor 3D PCAs
    if("Multiple_factor" %in% names(PCA_data) & grepl(pattern = 2,pca_type)){
      for(PCA_comp in names(PCA_data[["Multiple_factor"]])){
        colors_3d = brewer.pal(n = 9,name = "Set1")[match(PCA_data[["Multiple_factor"]][[PCA_comp]][["group"]],table = unique(PCA_data[["Multiple_factor"]][[PCA_comp]][["group"]]))]
        if(run_factor %in% attr(PCA_data[["Multiple_factor"]][[PCA_comp]], which = "factor") & !any(sapply(multi_done,function(x){all(x==attr(PCA_data[["Multiple_factor"]][[PCA_comp]], which = "factor"))}))){
          multi_done = append(multi_done,values = list(attr(PCA_data[["Multiple_factor"]][[PCA_comp]], which = "factor")))
          dir.create("animation_merge",showWarnings = FALSE)
          for(degree in 1:360) {
            open3d()
            par3d(windowRect = c(20, 30, 1080, 1080),dev = unname(rgl.dev.list()))
            plot3d(x = PCA_data[["Multiple_factor"]][[PCA_comp]][["PC1"]],
                   y = PCA_data[["Multiple_factor"]][[PCA_comp]][["PC2"]],
                   z = PCA_data[["Multiple_factor"]][[PCA_comp]][["PC3"]],
                   col = colors_3d,
                   size = 2,
                   xlab = paste0("PC1: ",round(100 * attr(PCA_data[["Multiple_factor"]][[PCA_comp]], "percentVar"))[1],"% variance"),
                   ylab = paste0("PC2: ",round(100 * attr(PCA_data[["Multiple_factor"]][[PCA_comp]], "percentVar"))[2],"% variance"),
                   zlab = paste0("PC3: ",round(100 * attr(PCA_data[["Multiple_factor"]][[PCA_comp]], "percentVar"))[3],"% variance"),
                   type = "s") +
              legend3d("topright", legend = unique(PCA_data[["Multiple_factor"]][[PCA_comp]][["group"]]), col = brewer.pal(n = 9,name = "Set1")[1:length(unique(PCA_data[["Multiple_factor"]][[PCA_comp]][["group"]]))], pch = 16, cex=1, inset=c(0.02))
            view3d(userMatrix=rotationMatrix(2*pi * degree/360, 0, 1, 0))
            rgl.snapshot(filename=paste("animation_merge/frame-",
                                        sprintf("%03d", degree), ".png", sep=""))
            while(length(rgl.dev.list()) != 0){rgl.close()}
            cat("\rRunning multiple factor ",PCA_comp," | ",format(round((degree/360)*100,digits = 2),nsmall = 2),"%", sep = "")
          }
          cat("\nRunning ffmpeg...\n")
          try(system(paste0("ffmpeg -y -hide_banner -loglevel warning -r 60 -y -i ",getwd(),"/animation_merge/frame-%03d.png ./",rlog_vst,"/PCA/PCA_3D_",rlog_vst,"_",genes_isoforms,"_",PCA_comp,".mp4")))
          unlink("animation_merge",recursive = TRUE)
        }
      }
    }
    
    # DEGs 3D PCAs
    if("DEGs" %in% names(PCA_data) & grepl(pattern = 3,pca_type)){
      if(run_factor %in% names(PCA_data[["DEGs"]])){
        for(PCA_comp in names(PCA_data[["DEGs"]][[run_factor]])){
          colors_3d = brewer.pal(n = 9,name = "Set1")[match(PCA_data[["DEGs"]][[run_factor]][[PCA_comp]][["group"]],table = unique(PCA_data[["DEGs"]][[run_factor]][[PCA_comp]][["group"]]))]
          dir.create("animation_merge",showWarnings = FALSE)
          for(degree in 1:360) {
            open3d()
            par3d(windowRect = c(20, 30, 1080, 1080),dev = unname(rgl.dev.list()))
            plot3d(x = PCA_data[["DEGs"]][[run_factor]][[PCA_comp]][["PC1"]],
                   y = PCA_data[["DEGs"]][[run_factor]][[PCA_comp]][["PC2"]],
                   z = PCA_data[["DEGs"]][[run_factor]][[PCA_comp]][["PC3"]],
                   col = colors_3d,
                   size = 2,
                   xlab = paste0("PC1: ",round(100 * attr(PCA_data[["DEGs"]][[run_factor]][[PCA_comp]], "percentVar"))[1],"% variance"),
                   ylab = paste0("PC2: ",round(100 * attr(PCA_data[["DEGs"]][[run_factor]][[PCA_comp]], "percentVar"))[2],"% variance"),
                   zlab = paste0("PC3: ",round(100 * attr(PCA_data[["DEGs"]][[run_factor]][[PCA_comp]], "percentVar"))[3],"% variance"),
                   type = "s") +
              legend3d("topright", legend = unique(PCA_data[["DEGs"]][[run_factor]][[PCA_comp]][["group"]]), col = brewer.pal(n = 9,name = "Set1")[1:length(unique(PCA_data[["DEGs"]][[run_factor]][[PCA_comp]][["group"]]))], pch = 16, cex=1, inset=c(0.02))
            view3d(userMatrix=rotationMatrix(2*pi * degree/360, 0, 1, 0))
            rgl.snapshot(filename=paste("animation_merge/frame-",
                                        sprintf("%03d", degree), ".png", sep=""))
            while(length(rgl.dev.list()) != 0){rgl.close()}
            cat("\rRunning DEGs ",PCA_comp," | ",format(round((degree/360)*100,digits = 2),nsmall = 2),"%", sep = "")
          }
          cat("\nRunning ffmpeg...\n")
          try(system(paste0("ffmpeg -y -hide_banner -loglevel warning -r 60 -y -i ",getwd(),"/animation_merge/frame-%03d.png ./",rlog_vst,"/PCA/PCA_3D_DEGs_",rlog_vst,"_",genes_isoforms,"_",PCA_comp,".mp4")))
          unlink("animation_merge",recursive = TRUE)
        }
      }
    }
  }
}

# Kmeans pheatmap function for seed setting
pheatmap_seed = function(mat,
                         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                         kmeans_k = NA,
                         breaks = NA,
                         border_color = "grey60",
                         cellwidth = NA,
                         cellheight = NA,
                         scale = "none",
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         clustering_method = "complete",
                         clustering_callback = identity2,
                         cutree_rows = NA,
                         cutree_cols = NA,
                         treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows, 50, 0), treeheight_col = ifelse((class(cluster_cols) == "hclust") || cluster_cols, 50, 0),
                         legend = TRUE,
                         legend_breaks = NA,
                         legend_labels = NA,
                         annotation_row = NA,
                         annotation_col = NA,
                         annotation = NA,
                         annotation_colors = NA,
                         annotation_legend = TRUE,
                         annotation_names_row = TRUE,
                         annotation_names_col = TRUE,
                         drop_levels = TRUE,
                         show_rownames = T,
                         show_colnames = T,
                         main = NA,
                         fontsize = 10,
                         fontsize_row = fontsize,
                         fontsize_col = fontsize,
                         angle_col = c("270", "0", "45", "90", "315"),
                         display_numbers = F,
                         number_format = "%.2f",
                         number_color = "grey30",
                         fontsize_number = 0.8 * fontsize,
                         gaps_row = NULL,
                         gaps_col = NULL,
                         labels_row = NULL,
                         labels_col = NULL,
                         filename = NA,
                         width = NA,
                         height = NA,
                         silent = FALSE,
                         na_col = "#DDDDDD",
                         seed_input = 123, ...){
  # Set labels
  if(is.null(labels_row)){
    labels_row = rownames(mat)
  }
  if(is.null(labels_col)){
    labels_col = colnames(mat)
  }
  # Preprocess matrix
  mat = as.matrix(mat)
  if(scale != "none"){
    mat = scale_mat(mat, scale)
    if(is.na2(breaks)){
      breaks = generate_breaks(mat, length(color), center = T)
    }
  }
  # Kmeans
  if(!is.na(kmeans_k)){
    # Cluster data
    set.seed(seed_input)
    km = kmeans(mat, kmeans_k, iter.max = 100)
    mat = km$centers
    # Compose rownames
    t = table(km$cluster)
    labels_row = sprintf("Cluster: %s Size: %d", names(t), t)
  }
  else{
    km = NA
  }
  # Format numbers to be displayed in cells
  if(is.matrix(display_numbers) | is.data.frame(display_numbers)){
    if(nrow(display_numbers) != nrow(mat) | ncol(display_numbers) != ncol(mat)){
      stop("If display_numbers provided as matrix, its dimensions have to match with mat")
    }
    display_numbers = as.matrix(display_numbers)
    fmat = matrix(as.character(display_numbers), nrow = nrow(display_numbers), ncol = ncol(display_numbers))
    fmat_draw = TRUE
  }else{
    if(display_numbers){
      fmat = matrix(sprintf(number_format, mat), nrow = nrow(mat), ncol = ncol(mat))
      fmat_draw = TRUE
    }else{
      fmat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
      fmat_draw = FALSE
    }
  }
  # Do clustering
  if((class(cluster_rows) == "hclust") || cluster_rows){
    if(class(cluster_rows) == "hclust"){
      tree_row = cluster_rows
    } else {
      tree_row = cluster_mat(mat, distance = clustering_distance_rows, method = clustering_method)
      tree_row = clustering_callback(tree_row, mat)
    }
    mat = mat[tree_row$order, , drop = FALSE]
    fmat = fmat[tree_row$order, , drop = FALSE]
    labels_row = labels_row[tree_row$order]
    if(!is.na(cutree_rows)){
      gaps_row = find_gaps(tree_row, cutree_rows)
    }else{
      gaps_row = NULL
    }
  }else{
    tree_row = NA
    treeheight_row = 0
  }
  if((class(cluster_cols) == "hclust") || cluster_cols){
    if(class(cluster_cols) == "hclust"){
      tree_col = cluster_cols
    } else {
      tree_col = cluster_mat(t(mat), distance = clustering_distance_cols, method = clustering_method)
      tree_col = clustering_callback(tree_col, t(mat))
    }
    mat = mat[, tree_col$order, drop = FALSE]
    fmat = fmat[, tree_col$order, drop = FALSE]
    labels_col = labels_col[tree_col$order]
    if(!is.na(cutree_cols)){
      gaps_col = find_gaps(tree_col, cutree_cols)
    }else{
      gaps_col = NULL
    }
  }else{
    tree_col = NA
    treeheight_col = 0
  }
  attr(fmat, "draw") = fmat_draw
  # Colors and scales
  if(!is.na2(legend_breaks) & !is.na2(legend_labels)){
    if(length(legend_breaks) != length(legend_labels)){
      stop("Lengths of legend_breaks and legend_labels must be the same")
    }
  }
  if(is.na2(breaks)){
    breaks = generate_breaks(as.vector(mat), length(color))
  }
  if(!is.infinite(min(breaks))){
    breaks = c(-Inf, breaks)
    color = c(color[1], color)
  }
  if(!is.infinite(max(breaks))){
    breaks = c(breaks, Inf)
    color = c(color, color[length(color)])
  }
  non_inf_breaks = breaks[!is.infinite(breaks)]
  if (legend & is.na2(legend_breaks)) {
    legend = grid.pretty(range(as.vector(non_inf_breaks)))
    names(legend) = legend
  }
  else if(legend & !is.na2(legend_breaks)){
    legend = legend_breaks[legend_breaks >= min(non_inf_breaks) & legend_breaks <= max(non_inf_breaks)]
    if(!is.na2(legend_labels)){
      legend_labels = legend_labels[legend_breaks >= min(non_inf_breaks) & legend_breaks <= max(non_inf_breaks)]
      names(legend) = legend_labels
    }
    else{
      names(legend) = legend
    }
  }
  else {
    legend = NA
  }
  mat = scale_colours(mat, col = color, breaks = breaks, na_col = na_col)
  # Preparing annotations
  if(is.na2(annotation_col) & !is.na2(annotation)){
    annotation_col = annotation
  }
  # Select only the ones present in the matrix
  if(!is.na2(annotation_col)){
    annotation_col = annotation_col[colnames(mat), , drop = F]
  }
  if(!is.na2(annotation_row)){
    annotation_row = annotation_row[rownames(mat), , drop = F]
  }
  annotation = c(annotation_row, annotation_col)
  annotation = annotation[unlist(lapply(annotation, function(x) !is.na2(x)))]
  if(length(annotation) != 0){
    annotation_colors = generate_annotation_colours(annotation, annotation_colors, drop = drop_levels)
  }
  else{
    annotation_colors = NA
  }
  if(!show_rownames){
    labels_row = NULL
  }
  if(!show_colnames){
    labels_col = NULL
  }
  # Set colname rotating parameters
  angle_col = as.character(angle_col)
  angle_col = match.arg(angle_col)
  if(angle_col == "0"){
    angle_col = 0
    hjust_col = 0.5
    vjust_col = 1
  }
  if(angle_col == "45"){
    angle_col = 45
    hjust_col = 1
    vjust_col = 1
  }
  if(angle_col == "90"){
    angle_col = 90
    hjust_col = 1
    vjust_col = 0.5
  }
  if(angle_col == "270"){
    angle_col = 270
    hjust_col = 0
    vjust_col = 0.5
  }
  if(angle_col == "315"){
    angle_col = 315
    hjust_col = 0
    vjust_col = 1
  }
  # Draw heatmap
  pdf(file = NULL)
  gt = heatmap_motor(mat,
                     border_color = border_color,
                     cellwidth = cellwidth,
                     cellheight = cellheight,
                     treeheight_col = treeheight_col,
                     treeheight_row = treeheight_row,
                     tree_col = tree_col,
                     tree_row = tree_row,
                     filename = filename,
                     width = width,
                     height = height,
                     breaks = breaks,
                     color = color,
                     legend = legend,
                     annotation_row = annotation_row,
                     annotation_col = annotation_col,
                     annotation_colors = annotation_colors,
                     annotation_legend = annotation_legend,
                     annotation_names_row = annotation_names_row,
                     annotation_names_col = annotation_names_col,
                     main = main,
                     fontsize = fontsize,
                     fontsize_row = fontsize_row,
                     fontsize_col = fontsize_col,
                     hjust_col = hjust_col,
                     vjust_col = vjust_col,
                     angle_col = angle_col,
                     fmat = fmat,
                     fontsize_number = fontsize_number,
                     number_color = number_color,
                     gaps_row = gaps_row,
                     gaps_col = gaps_col,
                     labels_row = labels_row,
                     labels_col = labels_col,
                     ...)
  dev.off()
  if(is.na(filename) & !silent){
    grid.newpage()
    grid.draw(gt)
  }
  invisible(structure(list(tree_row = tree_row, tree_col = tree_col, kmeans = km, gtable = gt), class = "pheatmap"))
}

# Set the environment for pheatmap_seed to that of pheatmap  
environment(pheatmap_seed) = environment(pheatmap)
}
#####    Start run    #####
  if("DESeq2" %in% section | "LRT" %in% section){
  cat("#####    Importing data    ######\n", sep = "")
    if(exists("input_path")){
      file_names = file.path(input_path,basename(file_names))
    }
  if(Mapper == "RSEM"){
  names(file_names) = rownames(experimental_design)
  if(genes_isoforms == "genes"){
    RSEM_counts = tximport(files = file_names,type = "rsem",txIn = FALSE,txOut = FALSE)
  }
  if(genes_isoforms == "isoforms"){
    RSEM_counts = tximport(files = file_names,type = "rsem",txIn = TRUE,txOut = TRUE)
  }
  RSEM_counts$length[RSEM_counts$length == 0] = 1
  data_set = DESeqDataSetFromTximport(txi = RSEM_counts,colData = experimental_design,design = design_formula)
  gene_length = RSEM_counts$length
  count_table = assay(data_set)
  }
    
    if(Mapper == "Kallisto"){
      names(file_names) = rownames(experimental_design)
      kallisto_counts = tximport(files = file_names,type = "kallisto", tx2gene = tx2gene)
      kallisto_counts$length[kallisto_counts$length == 0] = 1
      data_set = DESeqDataSetFromTximport(txi = kallisto_counts,colData = experimental_design,design = design_formula)
      gene_length = kallisto_counts$length
      count_table = assay(data_set)
    }
    
  if(Mapper == "HTseq-count"){
    sampleTable = data.frame(sampleName = rownames(experimental_design),
                             filename = basename(file_names))
    sampleTable = cbind(sampleTable,experimental_design)
    data_set = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = input_path,design = design_formula)
    count_table = assay(data_set)
    if(any(grepl(pattern = "^N_unmapped$|^N_multimapping$|^N_noFeature$|^N_ambiguous$",rownames(count_table)))){
      count_table = count_table[!grepl(pattern = "^N_unmapped$|^N_multimapping$|^N_noFeature$|^N_ambiguous$",rownames(count_table)),]
    }
  }
    
  if(Mapper == "Counts"){
    for(i in file_names){
      cat("Reading file: ",i,"\n", sep = "")
      if(i == file_names[1]){
        count_table = read.table(file = i,header = TRUE,sep = "\t")[c(count_genes,count_column)]
        rownames(count_table) = count_table[,count_genes]
        if(!exists("gtf")){
          gene_length = read.table(file = i,header = TRUE,sep = "\t")[count_length]
        if(length(count_column) > 1){
          gene_length[,2:length(count_column)] = gene_length[[1]]
        }
        rownames(gene_length) = count_table[,count_genes]
        }
        count_table = count_table[,-count_genes,drop = FALSE]
      }else{
        count_table = cbind(count_table,read.table(file = i,header = TRUE,sep = "\t")[count_column])
        if(!exists("gtf")){
          gene_length = cbind(gene_length,read.table(file = i,header = TRUE,sep = "\t")[count_length])
        }
      }
    }
        colnames(count_table) = rownames(experimental_design)
    data_set = DESeqDataSetFromMatrix(countData = count_table,
                                      colData = experimental_design,
                                      design = design_formula,
                                      tidy = FALSE,
                                      ignoreRank = FALSE)
  }
    if(NOISeq_correction){
      count_table = round(ARSyNseq(data = readData(data = filtered.data(dataset = count_table,
                                                                        depth = colSums(count_table),
                                                                        p.adj = "fdr",
                                                                        method = 3,
                                                                        norm = FALSE,
                                                                        factor = primary_factor),
                                                   factors = experimental_design,
                                                   length = as.matrix(gene_length)[,1]),
                                   factor = primary_factor,
                                   norm = "n")@assayData$exprs)
      gene_length = gene_length[rownames(count_table),,drop = FALSE]
      data_set = DESeqDataSetFromMatrix(countData = count_table,
                                        colData = experimental_design,
                                        design = design_formula,
                                        tidy = FALSE,
                                        ignoreRank = FALSE)
    }
      TPM = if(ncol(count_table) == ncol(gene_length)){
        as.data.frame(lapply(1:ncol(count_table),function(x){
          RPK = count_table[,x]/(gene_length[,x]/1000)
          return(RPK/(sum(RPK)/(10^6)))
      }),row.names = rownames(count_table),col.names = colnames(count_table))
      }else{
        as.data.frame(lapply(1:ncol(count_table),function(x){
          RPK = count_table[,x]/(gene_length[,1]/1000)
          return(RPK/(sum(RPK)/(10^6)))
        }),row.names = rownames(count_table),col.names = colnames(count_table))
      }
      write.table(TPM,file = "TPM.txt",quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
      write.table(TPM[rowSums(TPM) > 0,],file = "TPM_non_zero.txt",quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
      if("bcbioRNASeq" %in% rownames(installed.packages())){
        library(bcbioRNASeq)
        TMM = tmm(count_table)
        write.table(TMM,file = "TMM.txt",quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
        write.table(TMM[rowSums(TMM) > 0,],file = "TMM_non_zero.txt",quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
      }
    data_set = data_set[rowSums(counts(data_set)) >= 10]
    for(i in factors){
      data_set@colData[[i]] = relevel(data_set@colData[[i]], ref = variable_ref[i])
    }

  ##### Run DESeq2 #####
  cat("#####    Running DESeq2   ######\n")
  time_start=Sys.time()
  if("DESeq2" %in% section){
    data_set_DESeq = DESeq(data_set,parallel = TRUE)
    cat("DESeq ",genes_isoforms," completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }
  if("LRT" %in% section){
    data_set_DESeq = DESeq(data_set,test = "LRT",reduced = reduced_formula,parallel = TRUE)
    cat("LRT ",genes_isoforms," completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }
  
    save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_DESeq2_data.RData"))
  }
  
  if("Transformation" %in% section){
  #####   Transform counts   #####
  if(rlog_vst == "rlog"){
    cat("#####    rLog transforming ",genes_isoforms," counts    ######\n", sep = "")
    time_start=Sys.time()
    data_set_transform = rlog(data_set_DESeq,blind = FALSE)
    cat("rLog transforming ",genes_isoforms," counts completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }
  if(rlog_vst == "VST"){
    cat("#####    VST transforming ",genes_isoforms," counts    ######\n", sep = "")
    time_start=Sys.time()
    data_set_transform = varianceStabilizingTransformation(data_set_DESeq,blind = FALSE)
    cat("VST transforming ",genes_isoforms," counts completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }
    wgcna_data = if(rlog_vst == "rlog"){
      assay(varianceStabilizingTransformation(data_set_DESeq,blind = FALSE))
      }else{
        assay(data_set_transform)
      }
    saveRDS(object = experimental_design,file = paste0(Experiment_name,"_",genes_isoforms,"_","factors_for_WGCNA.rds"))
    saveRDS(object = wgcna_data,file = paste0(Experiment_name,"_",genes_isoforms,"_","VST_for_WGCNA.rds"))
    rm(wgcna_data)
  save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_","transformed_data.RData"))
  }
  
  if(any(c("Heatmaps","Dispersion estimates, PCAs and PCoAs","Optimize K-means","K-means clustering","Hierarchial clustering","MA-plot","DEG heatmap","goseq GO analysis","topGO analysis","Wordcloud","Venn diagram","Variable heatmap and report") %in% section)){
    #####   Data matrix   #####
    dir.create(rlog_vst,showWarnings = FALSE)
    dir.create(paste0(rlog_vst,"/Heatmaps"),showWarnings = FALSE)
    data_set_matrix = assay(data_set_transform)
    # if(nrow(data_set_matrix) > 65536){
    #   cat("Too many genes. Removing lowest ",nrow(data_set_matrix) - 65536," expressed genes\n",sep = "")
    #   data_set_matrix = data_set_matrix[names(head(sort(abs(rowSums(data_set_matrix)),decreasing = TRUE),n = 65536)),]
    # }
    data_set_matrix_scaled = scale(data_set_matrix, center=TRUE, scale=TRUE)
    save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_","analysis_data.RData"))
    
    if("Heatmaps" %in% section){
    cat("#####    Creating matrices and heatmaps    ######\n")
      time_start=Sys.time()
      pheatmap(
        if(nrow(data_set_matrix) > 65536){
        data_set_matrix_scaled[names(head(sort(abs(rowSums(data_set_matrix_scaled)),decreasing = TRUE),n = 65536)),]
      }else{
        data_set_matrix_scaled
      },
             show_rownames = FALSE,
             cluster_rows = heatmap_row_clust,
             annotation_col = as.data.frame(colData(data_set)),
             filename = paste0(rlog_vst,"/Heatmaps/Heatmap_",rlog_vst,"_",genes_isoforms,".png"),
             annotation_colors = color_select)
    if(exists("gene_subset")){
      gene_subset = read.table(file = gene_subset_path,sep = "\t",header = FALSE)
      # gene_subset$V2 = gsub(pattern = "[-]|[.]",replacement = "_",x = gene_subset$V2)
      # gene_subset$V2[grep(pattern = "L-arginine biosynthesis I",x = gene_subset$V2)] = "L-arginine biosynthesis I"
      # gene_subset$V2[grep(pattern = "glutaminyl-tRNA biosynthesis via transamidation",x = gene_subset$V2)] = "glutaminyl-tRNA biosynthesis"
      # gene_subset$V2[grep(pattern = "superpathway of branched chain amino acid",x = gene_subset$V2)] = "amino acid biosynthesis"
      rownames(gene_subset) = gene_subset$V1
      gene_subset = gene_subset[,-1,drop = FALSE]
      # colnames(gene_subset) = "Pathway"
      pheatmap(data_set_matrix_scaled[rownames(data_set_matrix_scaled) %in% rownames(gene_subset),],
               show_rownames = FALSE,
               annotation_row = gene_subset,
               annotation_col = as.data.frame(colData(data_set)),
               filename = paste0(rlog_vst,"/Heatmaps/Heatmap_gene_subset_",rlog_vst,"_",genes_isoforms,".png"),
               annotation_colors = color_select)
    }
    save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_","analysis_data.RData"))
    cat("Matrices and heatmaps completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
    }
    
  if("Dispersion estimates, PCAs and PCoAs" %in% section){
  #####   Dispersion estimates   #####
  time_start=Sys.time()
  cat("#####    Creating dispersion estimates    ######\n")
  png(filename = paste0("Plot_dispertion_estimates_",genes_isoforms,".png"),width = 1920,height = 1080,units = "px")
  plotDispEsts(data_set_DESeq)
  while (!is.null(dev.list())){dev.off()}
  cat("Dispersion estimates completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  #####   Single factor PCA   #####
  time_start=Sys.time()
  dir.create(paste0(rlog_vst,"/PCA"),showWarnings = FALSE)
  PCA_data = list()
  cat("#####    Creating PCA and factor heatmap plots    ######\n")
  for(i in factors){
      PCA_data[["Single_factor"]][[i]] = plotPCA_PC123(object = data_set_transform,intgroup=i,returnData = TRUE)
      if(!is.null(PCA_data[["Single_factor"]][[i]])){
      png(filename = paste0(rlog_vst,"/PCA/PCA_",rlog_vst,"_",genes_isoforms,"_",i,".png"),width = 1920,height = 1080,units = "px")
      plot_temp = ggplot(PCA_data[["Single_factor"]][[i]],
                     aes(PC1, PC2, color = eval(expr = parse(text = i)), group = experimental_design[[i]], label = PCA_data[["Single_factor"]][[i]][["name"]])) +
            geom_point(size=4) +
            xlab(paste0("PC1: ",round(100 * attr(PCA_data[["Single_factor"]][[i]], "percentVar"))[1],"% variance")) +
            ylab(paste0("PC2: ",round(100 * attr(PCA_data[["Single_factor"]][[i]], "percentVar"))[2],"% variance")) +
            theme(text = element_text(size = 20)) +
            coord_fixed() +
            scale_color_discrete(name = i) +
            geom_text_repel(size = 8,vjust = 0,nudge_y = 3,segment.size = 0, show.legend = FALSE,segment.color = "transparent") +
            guides(color=guide_legend(override.aes=list(fill=NA)))
      if(max(table(PCA_data[["Single_factor"]][[i]][["group"]])) > 3){
        plot_temp = plot_temp +
        stat_ellipse(geom = "polygon", alpha = 0.25, aes(fill = PCA_data[["Single_factor"]][[i]][["group"]]), lwd = 0, show.legend = any(table(PCA_data[["Single_factor"]][[i]][["group"]]) > 3)) +
        labs(fill = "Ellipse")
      }
      print(plot_temp)
      while (!is.null(dev.list())){dev.off()}
      
      
      png(filename = paste0(rlog_vst,"/PCA/PCA_",rlog_vst,"_",genes_isoforms,"_",i,"_PC2.png"),width = 1920,height = 1080,units = "px")
      plot_temp = ggplot(PCA_data[["Single_factor"]][[i]],
                     aes(PC2, PC3, color = eval(expr = parse(text = i)), group = experimental_design[[i]], label = PCA_data[["Single_factor"]][[i]][["name"]])) +
              geom_point(size=4) +
              xlab(paste0("PC2: ",round(100 * attr(PCA_data[["Single_factor"]][[i]], "percentVar"))[2],"% variance")) +
              ylab(paste0("PC3: ",round(100 * attr(PCA_data[["Single_factor"]][[i]], "percentVar"))[3],"% variance")) +
              theme(text = element_text(size = 20)) +
              coord_fixed() +
              scale_color_discrete(name = i) +
              geom_text_repel(size = 8,vjust = 0,nudge_y = 3,segment.size = 0, show.legend = FALSE,segment.color = "transparent") +
              guides(color=guide_legend(override.aes=list(fill=NA)))
      if(max(table(PCA_data[["Single_factor"]][[i]][["group"]])) > 3){
        plot_temp = plot_temp +
              stat_ellipse(geom = "polygon", alpha = 0.25, aes(fill = PCA_data[["Single_factor"]][[i]][["group"]]), lwd = 0, show.legend = any(table(PCA_data[["Single_factor"]][[i]][["group"]]) > 3)) +
              labs(fill = "Ellipse")
      }
      print(plot_temp)
      while (!is.null(dev.list())){dev.off()}
      }
  }
  
  #####   Multi-factor PCA and PCoA   #####
  if(length(factors) > 1){
  for(i in factors[-length(factors)]){
  for(j in factors[(which(factors %in% i)+1):length(factors)]){
    PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]] = plotPCA_PC123(object = data_set_transform,intgroup=c(i,j),returnData = TRUE)
    if(!is.null(PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]])){
    attr(PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]], which = "factor") = c(i,j)
    png(filename = paste0(rlog_vst,"/PCA/PCA_",rlog_vst,"_",genes_isoforms,"_",i,"_vs_",j,".png"),width = 1920,height = 1080,units = "px")
    plot_temp = ggplot(PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]],
                       aes(PC1, PC2, color = eval(expr = parse(text = i)), shape = eval(expr = parse(text = j)), label = PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]][["name"]])) +
            geom_point(size=4) +
            xlab(paste0("PC1: ",round(100 * attr(PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]], "percentVar"))[1],"% variance")) +
            ylab(paste0("PC2: ",round(100 * attr(PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]], "percentVar"))[2],"% variance")) +
            theme(text = element_text(size = 20)) +
            coord_fixed() +
            scale_shape_manual(values = 15:100,name = j) +
            scale_color_discrete(name = i) +
            geom_text_repel(size = 8,vjust = 0,nudge_y = 3,segment.size = 0, show.legend = FALSE,segment.color = "transparent") +
            guides(color=guide_legend(override.aes=list(fill=NA)))
      if(max(table(PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]][["group"]])) > 3){
        plot_temp = plot_temp +
        stat_ellipse(geom = "polygon",alpha = 0.5,aes(fill = PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]][["group"]]), lwd = 0, show.legend = any(table(PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]][["group"]]) > 3)) +
        labs(fill = "Ellipse")
      }
    print(plot_temp)
    while (!is.null(dev.list())){dev.off()}
    
    png(filename = paste0(rlog_vst,"/PCA/PCA_",rlog_vst,"_",genes_isoforms,"_",i,"_vs_",j,"_PC2.png"),width = 1920,height = 1080,units = "px")
    plot_temp = ggplot(PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]],
                   aes(PC2, PC3, color = eval(expr = parse(text = i)), shape = eval(expr = parse(text = j)), label = PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]][["name"]])) +
            geom_point(size=4) +
            xlab(paste0("PC2: ",round(100 * attr(PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]], "percentVar"))[2],"% variance")) +
            ylab(paste0("PC3: ",round(100 * attr(PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]], "percentVar"))[3],"% variance")) +
            theme(text = element_text(size = 20)) +
            coord_fixed() +
            scale_shape_manual(values = 15:100,name = j) +
            scale_color_discrete(name = i) +
            geom_text_repel(size = 8,vjust = 0,nudge_y = 3,segment.size = 0, show.legend = FALSE,segment.color = "transparent") +
            guides(color=guide_legend(override.aes=list(fill=NA)))
    if(max(table(PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]][["group"]])) > 3){
      plot_temp = plot_temp +      
      stat_ellipse(geom = "polygon", alpha = 0.5, aes(fill = PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]][["group"]]), lwd = 0, show.legend = any(table(PCA_data[["Multiple_factor"]][[paste0(i,"_vs_",j)]][["group"]]) > 3)) +
      labs(fill = "Ellipse")
    }
    print(plot_temp)
    while (!is.null(dev.list())){dev.off()}
    }
    
    data_set_deseq_pcoa = pcoa(vegdist(t(data_set_matrix),method="manhattan")/1000)
    pcoa_scores = data_set_deseq_pcoa$vectors
    png(filename = paste0(rlog_vst,"/PCA/PCoA_",rlog_vst,"_",genes_isoforms,"_",i,"_",j,".png"),width = 1920,height = 1080,units = "px")
    par(oma = c(4, 1, 1, 1))
    plot(pcoa_scores[,"Axis.1"],
         pcoa_scores[,"Axis.2"],
         xlab = "PCo1",
         ylab = "PCo2",
         main = paste0("Principal Coordinate Analysis of ",i," and ",j),
         type = "p",
         pch = match(x = experimental_design[[i]],table = unique(experimental_design[[i]])),
         col = match(x = experimental_design[[j]],table = unique(experimental_design[[j]])),
         cex = 1.5,
         cex.axis = 1.5,
         cex.lab = 1.5,
         cex.main = 1.5,
         cex.sub = 1.5)
    text(pcoa_scores[,"Axis.1"],
         pcoa_scores[,"Axis.2"],
         rownames(data_set_deseq_pcoa$vectors),
         cex=1.5,
         pos=4)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom",
           inset = c(0,0),
           cex = 1.5,
           bty = 'n',
           horiz = TRUE,
           xpd = TRUE,
           legend = paste0(rep(unique(experimental_design[[i]]), each = length(unique(experimental_design[[j]])))," / ",rep(unique(experimental_design[[j]]), times = length(unique(experimental_design[[i]])))),
           pch = rep(seq_along(unique(experimental_design[[i]])),each = length(unique(experimental_design[[j]]))),
           col = rep(seq_along(unique(experimental_design[[j]])),length(unique(experimental_design[[i]]))))
    while (!is.null(dev.list())){dev.off()}
  }
  }
  }
  save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_","analysis_data.RData"))
  cat("PCAs and PCoAs completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }
    
  #######      Calculate optimal K      #######
  if("Optimize K-means" %in% section){
    time_start=Sys.time()
    cat("#####    Calculating optimal K   ######\n")
    cpucores = threads
    options("mc.cores" = cpucores)
    registerDoParallel(cpucores)
    k_clusGap = clusGapKB(x = data_set_matrix_scaled,FUNcluster = kmeans,iter.max = 50,K.max = 50,B = if(exists("k_boot_max")){k_boot_max}else{100})
    k_clusters = with(k_clusGap, maxSE(Tab[,"gap"], Tab[,"SE.sim"], method="firstSEmax"))
    png(filename = paste0("./clusGap_k_",k_clusters,"_",rlog_vst,"_",genes_isoforms,".png"),width = 1920,height = 1080,units = "px")
    plot(k_clusGap, frame = FALSE, xlab = "k clusters", ylab = "Gap statistic",main = paste0("k clusters vs. gap (k = ",k_clusters,")"))
    abline(v = k_clusters)
    while (!is.null(dev.list())){dev.off()}

    cat("Using K = ",k_clusters,"\n", sep = "")
    cat("Completed K optimization in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
    save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_","analysis_data.RData"))
  }
  
  #######      K-means clustering      #########
  if("K-means clustering" %in% section){
    cat("#####    K-means clustering and creating matplots    ######\n")
    time_start=Sys.time()
      set.seed(random_seed)
      km.res = kmeans(x = data_set_matrix_scaled,centers = k_clusters, iter.max = 100)
      samples_df = setDT(as.data.frame(data_set_matrix_scaled), keep.rownames = TRUE)[]
      km.table = setDT(as.data.frame(km.res$cluster), keep.rownames = TRUE)[]
      colnames(km.table)[2]="Cluster"
      write.table(km.table,file = "kmeans_cluster_table.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
      for(i in 1:k_clusters){
        cluster_exp = samples_df[samples_df$rn %in% km.table[km.table$Cluster==i,]$rn,]
        cluster_exp = cluster_exp[,-1]
        png(filename = paste0(rlog_vst,"/Matplot_Cluster_",i,"_",rlog_vst,"_",genes_isoforms,".png"),width = 1920,height = 1080,units = "px")
        par(mfrow=c(1,1))
        matplot(t(cluster_exp),
                type = "l",
                pch = 1,
                ylim = c(floor(min(cluster_exp[,2:ncol(cluster_exp)]))+floor((min(cluster_exp[,2:ncol(cluster_exp)])%%1)*100)/100,floor(max(cluster_exp[,2:ncol(cluster_exp)]))+ceiling((max(cluster_exp[,2:ncol(cluster_exp)])%%1)*100)/100),
                xaxt = "n",
                main = paste0("Cluster ",i," (",nrow(cluster_exp),")"),
                ylab = "z-score",
                col = gray.colors(12))
        lines(colMeans(cluster_exp, na.rm = FALSE, dims = 1), col="black", lwd=3)
        axis(1, at = 1:(ncol(cluster_exp)), labels = colnames(cluster_exp), cex.axis = 1)
        while (!is.null(dev.list())){dev.off()}
      }
        cat("#####    Creating heatmap for K clusters     ######\n")
      pheatmap_seed(data_set_matrix_scaled,
                    show_rownames = TRUE,
                    annotation_col = as.data.frame(colData(data_set)),
                    filename = paste0(rlog_vst,"/Heatmaps/Heatmap_kclusters","_",k_clusters,"_",rlog_vst,"_",genes_isoforms,".png"),
                    kmeans_k = k_clusters,
                    annotation_colors = color_select,
                    seed_input = random_seed)
      save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_","analysis_data.RData"))
      cat("Completed K-means clustering in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }
  
  ##### Create and plot Hierarchical clusters #####
  if("Hierarchial clustering" %in% section){
    cat("#####    Hierarchial clustering    ######\n")
    time_start=Sys.time()
    hclusters = hclust(dist(
      if(nrow(data_set_matrix) > 65536){
        data_set_matrix_scaled[names(head(sort(abs(rowSums(data_set_matrix_scaled)),decreasing = TRUE),n = 65536)),]
      }else{
        data_set_matrix
      }, method = "euclidean"), method = "complete")
    png(filename = paste0("./Hierarchical_Clustering","_",rlog_vst,"_",genes_isoforms,".png"),width = 1920,height = 1080,units = "px")
    plot(hclusters, labels = FALSE, hang = -1,ann = FALSE)
    title(xlab = "Genes",ylab = "Euclidean distance")
    if(exists("k_clusters")){
      rh = rect.hclust(hclusters, k = k_clusters, border = 1:k_clusters,cluster = cutree(hclusters, k = k_clusters))
      rh_lengths = head(cumsum(c(1, lengths(rh))), -1)
      rh_y = c(rep(-1,times = k_clusters))
      for(i in 2:k_clusters){
        if(rh_lengths[[i-1]] > rh_lengths[[i]] - 50){
          rh_y[[i]] = rh_y[[i-1]] + 0.5
          }
      }
      text(x = rh_lengths+35, y = rh_y, col="black", labels=1:k_clusters, font=1)
    }
    while (!is.null(dev.list())){dev.off()}
    cat("Completed Hierarchical clustering in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
    save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_","analysis_data.RData"))
  }
  
  #####   MA summary, MA plots, Gene heat map and GO annotation #####
  if("MA-plot" %in% section){
    time_start=Sys.time()
    for(i in resultsNames(data_set_DESeq)[-1]){
    cat("#####    Shrinking results and calculating ,MA summary and plots    ######\n")
    deseq_results_LFC = lfcShrink(data_set_DESeq,coef = i,type = "apeglm",parallel = TRUE)
    png(filename = paste0("./MA-plot_",i,"_LFC_deseq_results_",genes_isoforms,"_padj.png"),width = 1920,height = 1080,units = "px")
    plotMA(deseq_results_LFC, alpha = alpha, main = paste0(i), ylim = c(min(deseq_results_LFC$log2FoldChange),max(deseq_results_LFC$log2FoldChange)))
    while (!is.null(dev.list())){dev.off()}
    }
    save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_","analysis_data.RData"))
    cat("MA summary and MA plots completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }
    
  # if(remove_isoforms){
  #   gene_length = matrix(data = unlist(bplapply(X = unique(sub(pattern = isoform_pattern,replacement = "",x = rownames(gene_length))),FUN = function(x, y = gene_length){
  #     return(colMaxs(gene_length[grepl(pattern = x,x = rownames(gene_length)),,drop = FALSE]))
  #   })),ncol = ncol(gene_length),dimnames = list(unique(sub(pattern = isoform_pattern,replacement = "",x = rownames(gene_length))),colnames(gene_length)))
  #   save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_","analysis_data.RData"))
  # }
    
  deseq_results = list()
  deseq_results_sig = list()

  ## Individual factors
  for(f in factors){
    control_var = levels(data_set_DESeq@colData[[f]])[1]
    group_var = levels(data_set_DESeq@colData[[f]])[-1]
    for(i in group_var){
      time_start=Sys.time()
      cat("#####    Getting results for Factor: ",f," / Group: ",i,"   ######\n",sep = "")
      compare_var = paste0(f,"_",i,"_vs_",control_var)
      temp_results = results(data_set_DESeq,
                             contrast = c(f,i,control_var),
                             parallel = TRUE)
      sink(file = paste0(genes_isoforms,"_summary.txt"),append = TRUE)
      cat(i," vs. ",control_var, sep = "")
      summary(object = temp_results,alpha = alpha)
      sink()
      deseq_results[[compare_var]] = as.data.frame(temp_results)
      attr(deseq_results[[compare_var]],which = "factor") = f
    }
    rm(temp_results)
  }
  
  ## Interactions
  if(length(grep(pattern = ":",design_formula))){
  for(i in grep(pattern = "_vs_",resultsNames(data_set_DESeq)[-1],value = TRUE,invert = TRUE)){
    cat("#####    Getting results for interaction: ",i,"   ######\n",sep = "")
    temp_results = results(data_set_DESeq, name = i,parallel = TRUE)
    sink(file = paste0(genes_isoforms,"_summary.txt"),append = TRUE)
    cat(i,"_interaction", sep = "")
    summary(object = temp_results,alpha = alpha)
    sink()
    deseq_results[[i]] = as.data.frame(temp_results)
    attr(deseq_results[[compare_var]],which = "factor") = i
    rm(temp_results)
  }
  }
  save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_","analysis_data.RData"))
  
  dir.create(paste0(rlog_vst,"/Significant_DEGs"),showWarnings = FALSE)
  dir.create(paste0(rlog_vst,"/PCA"),showWarnings = FALSE)
  if(!exists("PCA_data")){
    PCA_data = list()
  }
  for(compare_var in names(deseq_results)){
      deseq_results[[compare_var]] = deseq_results[[compare_var]][order(deseq_results[[compare_var]][["padj"]]),]
      cat("Note: Removing ",sum(is.na(deseq_results[[compare_var]]))," empty genes from analysis\n", sep = "")
      deseq_results[[compare_var]] = deseq_results[[compare_var]][!is.na(deseq_results[[compare_var]][["padj"]]),]
      write.table(deseq_results[[compare_var]], file = paste0(compare_var,"_deseq_results_",genes_isoforms,"_padj.txt"),quote = FALSE,sep = "\t")
      
      #####    Create count matrix heat map    ######
      deseq_results_sig[[compare_var]] = deseq_results[[compare_var]][deseq_results[[compare_var]][["padj"]] < alpha & abs(deseq_results[[compare_var]][["log2FoldChange"]]) >= 1,]
      deseq_sig = rownames(deseq_results_sig[[compare_var]])
      DE_genes_sig = deseq_sig
      if(remove_isoforms){
        DE_genes_sig = unique(sub(pattern = isoform_pattern,replacement = "",x = DE_genes_sig))
      }
      write.table(x = DE_genes_sig,
                  file = paste0(rlog_vst,"/Significant_DEGs/Significant_DEGs_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_all.txt"),
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = FALSE)
      if(exists("k_clusters")){
        for(k in 1:k_clusters){
        DE_genes_sig = deseq_sig[deseq_sig %in% km.table[km.table$Cluster == k,]$rn]
        if(remove_isoforms){
          DE_genes_sig = unique(sub(pattern = isoform_pattern,replacement = "",x = DE_genes_sig))
        }
        write.table(x = DE_genes_sig,
                    file = paste0(rlog_vst,"/Significant_DEGs/Significant_DEGs_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,".txt"),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE)
        }
      }
      if("DEG heatmap" %in% section & length(deseq_sig) > 0){
      cat("#####    Creating DEG count matrix heat map    ######\n")
      pheatmap(data_set_matrix_scaled[deseq_sig,,drop = FALSE],
               show_rownames = FALSE,
               cluster_rows = heatmap_row_clust,
               annotation_col = as.data.frame(colData(data_set)),
               filename = paste0(rlog_vst,"/Heatmaps/Heatmap_significant_DEGs_",compare_var,"_",rlog_vst,"_",genes_isoforms,".png"),
               annotation_colors = color_select)
      if(exists("k_clusters")){
        if(length(deseq_sig) >= k_clusters){
        pheatmap_seed(data_set_matrix_scaled[deseq_sig,,drop = FALSE],
                      show_rownames = TRUE,
                      annotation_col = as.data.frame(colData(data_set)),
                      filename = paste0(rlog_vst,"/Heatmaps/Heatmap_significant_DEGs_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_k",k_clusters,".png"),
                      kmeans_k = k_clusters,
                      annotation_colors = color_select,
                      seed_input = random_seed)
        }}
        save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_","analysis_data.RData"))
        cat(" DEG heatmap completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
      }
      
      ######       PCA DEGs       ######
      if("Dispersion estimates, PCAs and PCoAs" %in% section){
      if(!is.null(attr(deseq_results[[compare_var]],which = "factor")) & length(deseq_sig) > 0){
      PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]] = plotPCA_PC123(object = data_set_transform[rownames(deseq_results_sig[[compare_var]]),],intgroup=attr(deseq_results[[compare_var]],which = "factor"),returnData = TRUE)
      if(!is.null(PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]])){
      png(filename = paste0(rlog_vst,"/PCA/PCA_DEGs_",rlog_vst,"_",genes_isoforms,"_",compare_var,".png"),width = 1920,height = 1080,units = "px")
      plot_temp = ggplot(PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]],
                         aes(PC1, PC2, color = eval(expr = parse(text = attr(deseq_results[[compare_var]],which = "factor"))), group = experimental_design[[attr(deseq_results[[compare_var]],which = "factor")]], label = PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]][["name"]])) +
        geom_point(size=4) +
        xlab(paste0("PC1: ",round(100 * attr(PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]], "percentVar"))[1],"% variance")) +
        ylab(paste0("PC2: ",round(100 * attr(PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]], "percentVar"))[2],"% variance")) +
        theme(text = element_text(size = 20)) +
        coord_fixed() +
        scale_color_discrete(name = attr(deseq_results[[compare_var]],which = "factor")) +
        geom_line(size = 0) +
        geom_text_repel(size = 8,vjust = 0,nudge_y = 3,segment.size = 0, show.legend = FALSE) +
        guides(color=guide_legend(override.aes=list(fill=NA)))
      if(max(table(PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]][["group"]])) > 3){
        plot_temp = plot_temp +
          stat_ellipse(geom = "polygon", alpha = 0.25, aes(fill = PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]][["group"]]), lwd = 0, show.legend = any(table(PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]][["group"]]) > 3)) +
          labs(fill = "Ellipse")
      }
      print(plot_temp)
      while (!is.null(dev.list())){dev.off()}
    
      png(filename = paste0(rlog_vst,"/PCA/PCA_DEGs_",rlog_vst,"_",genes_isoforms,"_",compare_var,"_PC2.png"),width = 1920,height = 1080,units = "px")
      plot_temp = ggplot(PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]],
                         aes(PC2, PC3, color = eval(expr = parse(text = attr(deseq_results[[compare_var]],which = "factor"))), group = experimental_design[[attr(deseq_results[[compare_var]],which = "factor")]], label = PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]][["name"]])) +
        geom_point(size=4) +
        xlab(paste0("PC2: ",round(100 * attr(PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]], "percentVar"))[2],"% variance")) +
        ylab(paste0("PC3: ",round(100 * attr(PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]], "percentVar"))[3],"% variance")) +
        theme(text = element_text(size = 20)) +
        coord_fixed() +
        scale_color_discrete(name = attr(deseq_results[[compare_var]],which = "factor")) +
        geom_text_repel(size = 8,vjust = 0,nudge_y = 3,segment.size = 0, show.legend = FALSE) +
        guides(color=guide_legend(override.aes=list(fill=NA)))
      if(max(table(PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]][["group"]])) > 3){
        plot_temp = plot_temp +
          stat_ellipse(geom = "polygon", alpha = 0.25, aes(fill = PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]][["group"]]), lwd = 0, show.legend = any(table(PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]][["group"]]) > 3)) +
          labs(fill = "Ellipse")
      }
      print(plot_temp)
      while (!is.null(dev.list())){dev.off()}
      }
      }
      }
      
      #######     goseq GO annotation      #########
      if("goseq GO analysis" %in% section){
        time_start=Sys.time()
        cat("#####    Starting goseq GO annotation    ######\n")
        dir.create(paste0(rlog_vst,"/goseq_GO_annotation/Annotation_table_sig_info"),showWarnings = FALSE, recursive = TRUE)
        dir.create(paste0(rlog_vst,"/goseq_GO_annotation/Annotation_table_sig"),showWarnings = FALSE)
        dir.create(paste0(rlog_vst,"/goseq_GO_annotation/Bootstrap"),showWarnings = FALSE)
        dir.create(paste0(rlog_vst,"/goseq_GO_annotation/Enrichment"),showWarnings = FALSE)
        dir.create(paste0(rlog_vst,"/goseq_GO_annotation/PWF"),showWarnings = FALSE)
        dir.create(paste0(rlog_vst,"/goseq_GO_annotation/Wordcloud"),showWarnings = FALSE)
        for(k in if(exists("k_clusters")){
                  c("all",1:k_clusters)
                  }else{
                  "all"}){
          if(k == "all"){
            DE_genes_sig = deseq_sig
            if(remove_isoforms){
              DE_genes_sig = unique(sub(pattern = isoform_pattern,replacement = "",x = DE_genes_sig))
            }
            write.table(x = GO_table[GO_table[["Gene"]] %in% DE_genes_sig,],
                        file = paste0(rlog_vst,"/goseq_GO_annotation/Annotation_table_sig_info/Annotation_table_sig_info_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_all.txt"),
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        sep = "\t")
            write.table(x = GO_table[GO_table[["Gene"]] %in% DE_genes_sig,c("Gene","GO_ID")],
                        file = paste0(rlog_vst,"/goseq_GO_annotation/Annotation_table_sig/Annotation_table_sig_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_all.txt"),
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        sep = "\t")
          }else{
            DE_genes_sig = deseq_sig[deseq_sig %in% km.table[km.table$Cluster == k,]$rn]
            if(remove_isoforms){
              DE_genes_sig = unique(sub(pattern = isoform_pattern,replacement = "",x = DE_genes_sig))
            }
            write.table(x = GO_table[GO_table[["Gene"]] %in% DE_genes_sig,],
                        file = paste0(rlog_vst,"/goseq_GO_annotation/Annotation_table_sig_info/Annotation_table_sig_info_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,".txt"),
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        sep = "\t")
            write.table(x = GO_table[GO_table[["Gene"]] %in% DE_genes_sig,c("Gene","GO_ID")],
                        file = paste0(rlog_vst,"/goseq_GO_annotation/Annotation_table_sig/Annotation_table_sig_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,".txt"),
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        sep = "\t")
          }
        
        all_genes = rownames(deseq_results[[compare_var]])
        if(remove_isoforms){
          rownames(gene_length) = sub(pattern = isoform_pattern,replacement = "",x = rownames(gene_length))
          all_genes = unique(sub(pattern = isoform_pattern,replacement = "",x = all_genes))
        }
        gene_vector = as.integer(all_genes %in% DE_genes_sig)
        
        # Create vector of gene lengths
        names(gene_vector) = all_genes
        gene_length_vector = rowMeans(gene_length[match(all_genes,rownames(gene_length)),,drop = FALSE])

        # Create Probability Weighting Function (PWF)
        png(filename = paste0(rlog_vst,"/goseq_GO_annotation/PWF/PWF_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,".png"),width = 1920,height = 1080,units = "px")
        pwf = nullp(DEgenes = gene_vector,bias.data = gene_length_vector)
        while (!is.null(dev.list())){dev.off()}
        rownames(pwf) = names(gene_length_vector)
        write.table(x = pwf,file = paste0(rlog_vst,"/goseq_GO_annotation/PWF/PWF_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,".txt"),quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
        
        # Estimate p-values for length-corrected DE genes
        pval = goseq(pwf = pwf,gene2cat = GO_table,method = "Wallenius")
        write.table(x = pval,file = paste0(rlog_vst,"/goseq_GO_annotation/Enrichment/GO_enriched_table_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,".txt"),quote = FALSE,row.names = FALSE)
        
        # Bootstrap for comparison to p-value estimation
        GO_sample = goseq(pwf = pwf,gene2cat = GO_table,method = "Sampling",repcnt = 2000)
        png(filename = paste0(rlog_vst,"/goseq_GO_annotation/Bootstrap/GO_bootstrap_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,".png"),width = 1920,height = 1080,units = "px")
        plot(log10(pval[,2]), log10(GO_sample[match(pval[,1],GO_sample[,1]),2]),xlab=expression('Log'[10]*' (Wallenius p-values)'),ylab=expression('Log'[10]*' (Sampling p-values)'),xlim=c(-3,0))
        abline(0,1,col = 3,lty = 2)
        while (!is.null(dev.list())){dev.off()}
        
        # Export table of FDR-corrected GO annotations
        write.table(x = pval$category[p.adjust(pval$over_represented_pvalue,method="BH") < FDR_cutoff],
                    file = paste0(rlog_vst,"/goseq_GO_annotation/Enrichment/GO_over_enriched_category_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,".txt"),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE)
        write.table(x = pval$term[p.adjust(pval$over_represented_pvalue,method="BH") < FDR_cutoff],
                    file = paste0(rlog_vst,"/goseq_GO_annotation/Enrichment/GO_over_enriched_term_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,".txt"),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE)
        write.table(x = pval$category[p.adjust(pval$under_represented_pvalue,method="BH") < FDR_cutoff],
                    file = paste0(rlog_vst,"/goseq_GO_annotation/Enrichment/GO_under_enriched_category_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,".txt"),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE)
        write.table(x = pval$term[p.adjust(pval$under_represented_pvalue,method="BH") < FDR_cutoff],
                    file = paste0(rlog_vst,"/goseq_GO_annotation/Enrichment/GO_under_enriched_term_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,".txt"),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE)
        
        }
        cat("goseq GO analysis completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
      }
      
      #######     topGO annotation      #########
      if("topGO analysis" %in% section){
        dir.create(paste0(rlog_vst,"/top_GO_annotation/topGO_graphs"),showWarnings = FALSE, recursive = TRUE)
        dir.create(paste0(rlog_vst,"/top_GO_annotation/topGO_DEGs"),showWarnings = FALSE)
        time_start=Sys.time()
        cat("#####    Starting topGO annotation    ######\n")
        for(k in if(exists("k_clusters")){
          c("all",1:k_clusters)
        }else{
          "all"}){
          if(k == "all"){
            DE_genes_sig = deseq_sig
            if(remove_isoforms){
              DE_genes_sig = unique(sub(pattern = isoform_pattern,replacement = "",x = DE_genes_sig))
            }
          }else{
            DE_genes_sig = deseq_sig[deseq_sig %in% km.table[km.table$Cluster == k,]$rn]
            if(remove_isoforms){
              DE_genes_sig = unique(sub(pattern = isoform_pattern,replacement = "",x = DE_genes_sig))
            }
          }
          if(sum(names(geneGO) %in% DE_genes_sig) > 0){
            geneList = factor(as.integer(names(geneGO) %in% DE_genes_sig),levels = c(0,1))
            names(geneList) = names(geneGO)
            for(ontology in c("BP","MF","CC")){
              sampleGOdata = new("topGOdata",
                                 description = "Simple session",
                                 ontology = ontology,
                                 allGenes = geneList,
                                 annot = annFUN.gene2GO,
                                 gene2GO = geneGO)
              enrich_result = runTest(sampleGOdata, statistic = "fisher")
              printGraph(object = sampleGOdata,
                         result = enrich_result,
                         firstSigNodes = 10,
                         fn.prefix = paste0(rlog_vst,"/top_GO_annotation/topGO_graphs/",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,"_",ontology),
                         useInfo = "all",
                         pdfSW = FALSE)
              
              if(sum(enrich_result@score <= alpha) > 0){
              GO_DEGs_df = as.data.frame(matrix(unlist(lapply(names(enrich_result@score[enrich_result@score <= alpha]),function(x){
                return(c(attr(GO_terms[[x]],which = "Term",exact = TRUE),
                         attr(GO_terms[[x]],which = "Definition",exact = TRUE)))
              })),ncol = 2,byrow = TRUE),row.names = names(names(enrich_result@score[enrich_result@score <= alpha])))
              colnames(GO_DEGs_df) = c("Term","Definition")
              write.table(GO_DEGs_df,file = paste0(rlog_vst,"/top_GO_annotation/topGO_DEGs/",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,"_",ontology,"_sig.txt"),quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
              
              if("Wordcloud" %in% section){
              cat("Creating enriched DEG Wordclouds...\n",sep = "")
              DEG_table = table(GO_DEGs_df[,"Term"])
              DEG_table = DEG_table[grep(pattern = "biological_process|cellular_component|molecular_function|*",x = names(DEG_table),invert = TRUE)]
              cat(format(round((3*match(k,if(exists("k_clusters")){c("all",1:k_clusters)}else{"all"}))/(3*if(exists("k_clusters")){k_clusters + 1}else{1}),digits = 2)*100,nsmall = 0),"% --> Comparison: ",compare_var," | k: ",k," | Ontology: ",ontology,"           \r",sep = "")
              if(length(DEG_table) > 0){
                png(filename = paste0(rlog_vst,"/Wordcloud/Enriched_Wordcloud_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,"_",ontology,".png"),width = 1080,height = 1080,units = "px")
                tryCatch(expr = {
                  wordcloud(words = names(DEG_table),
                            freq = DEG_table,
                            min.freq=1,
                            max.words=100,
                            random.order=FALSE,
                            rot.per=0.35,
                            use.r.layout=FALSE,
                            colors=brewer.pal(8, "Dark2"),
                            scale=c(5,0.4)
                  )
                },error = function(e) e)
                while (!is.null(dev.list())){dev.off()}
              }
              }
              }
            }
          }
        }
        # try(system(command = paste0("for i in ",rlog_vst,"/top_GO_annotation/topGO_graphs/*.ps; do convert -density 800 -rotate 90 $i ${i/.ps/.png} && rm -f $i; done")))
        cat("topGO analysis completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
      }
      
      ##### Create Wordcloud #####
      if("Wordcloud" %in% section){
        dir.create(paste0(rlog_vst,"/Wordcloud"),showWarnings = FALSE)
        cat("Creating DEG Wordclouds...\n",sep = "")
        for(k in if(exists("k_clusters")){
          c("all",1:k_clusters)
        }else{
          "all"}){
          if(k == "all"){
            DE_genes_sig = deseq_sig
            if(remove_isoforms){
              DE_genes_sig = unique(sub(pattern = isoform_pattern,replacement = "",x = DE_genes_sig))
            }
          }else{
            DE_genes_sig = deseq_sig[deseq_sig %in% km.table[km.table$Cluster == k,]$rn]
            if(remove_isoforms){
              DE_genes_sig = unique(sub(pattern = isoform_pattern,replacement = "",x = DE_genes_sig))
            }
          }
        for(ontology in c("BP","MF","CC")){
          DEG_table = GO_table[GO_table[["Gene"]] %in% DE_genes_sig,]
          DEG_table = DEG_table[!is.na(DEG_table[["Term"]]),,drop = FALSE]
          DEG_table = table(DEG_table[DEG_table[["Ontology"]] == ontology,"Term"])
          DEG_table = DEG_table[grep(pattern = "biological_process|cellular_component|molecular_function",x = names(DEG_table),invert = TRUE)]
          cat(format(round((3*match(k,if(exists("k_clusters")){c("all",1:k_clusters)}else{"all"}))/(3*if(exists("k_clusters")){k_clusters + 1}else{1}),digits = 2)*100,nsmall = 0),"% --> Comparison: ",compare_var," | k: ",k," | Ontology: ",ontology,"           \r",sep = "")
          if(length(DEG_table) > 0){
            png(filename = paste0(rlog_vst,"/Wordcloud/Wordcloud_",compare_var,"_",rlog_vst,"_",genes_isoforms,"_kcluster_",k,"_",ontology,".png"),width = 1080,height = 1080,units = "px")
            tryCatch(expr = {
              wordcloud(words = names(DEG_table),
                        freq = DEG_table,
                        min.freq=1,
                        max.words=100,
                        random.order=FALSE,
                        rot.per=0.35,
                        use.r.layout=FALSE,
                        colors=brewer.pal(8, "Dark2"),
                        scale=c(5,0.4)
              )
            },error = function(e) e)
            while (!is.null(dev.list())){dev.off()}
          }
        }
        }
        cat("\n")
      }
  save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_","analysis_data.RData"))
  }
  
    
  if("Venn diagram" %in% section){
  ##### Venn diagram #####
    # Intersect function
    dir.create(paste0(rlog_vst,"/Venn_Significant_DEGs"),showWarnings = FALSE)
    venn_calculate = function(x = NULL){
      intersect_list = list()
      intersect_list[[1]] = x
      for(i in 1:length(intersect_list[[1]])){
        attr(x = intersect_list[[1]][[i]],which = "Sample") = names(intersect_list[[1]][i])
      }
      for(intersect_order in 2:length(intersect_list[[1]])){
        intersect_list[[intersect_order]] = list()
        new_sample = 0
        for(sample_query in 1:length(intersect_list[[intersect_order-1]])){
          if(max(which(x = names(intersect_list[[1]]) %in% attr(intersect_list[[intersect_order-1]][[sample_query]],which = "Sample"))) < length(intersect_list[[1]])){
            for(sample_subject in (max(which(x = names(intersect_list[[1]]) %in% attr(intersect_list[[intersect_order-1]][[sample_query]],which = "Sample")))+1):length(intersect_list[[1]])){
              if(!attr(intersect_list[[1]][[sample_subject]],which = "Sample") %in% attr(intersect_list[[intersect_order-1]][[sample_query]],which = "Sample")){
                new_sample = new_sample+1
                intersect_list[[intersect_order]][[new_sample]] = intersect(
                  x = intersect_list[[intersect_order-1]][[sample_query]],
                  y = intersect_list[[1]][[sample_subject]]
                )
                attr(x = intersect_list[[intersect_order]][[new_sample]],which = "Sample") = c(
                  attr(intersect_list[[intersect_order-1]][[sample_query]],which = "Sample"),
                  attr(intersect_list[[1]][[sample_subject]],which = "Sample")
                )
              }
            }
          }
        }
      }
      venn_list = list()
      for(intersect_order in 1:(length(intersect_list)-1)){
        venn_list[[intersect_order]] = intersect_list[[intersect_order]]
        for(sample_query in 1:length(intersect_list[[intersect_order]])){
          for(sample_subject in 1:length(intersect_list[[intersect_order+1]])){
            if(any(attr(intersect_list[[intersect_order]][[sample_query]],which = "Sample") %in% attr(intersect_list[[intersect_order+1]][[sample_subject]],which = "Sample"))){
              venn_list[[intersect_order]][[sample_query]] = venn_list[[intersect_order]][[sample_query]][!venn_list[[intersect_order]][[sample_query]] %in% intersect_list[[intersect_order+1]][[sample_subject]]]
              attr(venn_list[[intersect_order]][[sample_query]], which = "Sample") = attr(intersect_list[[intersect_order]][[sample_query]],which = "Sample")
            }
          }
        }
      }
      venn_list[[length(intersect_list)]] = intersect_list[[length(intersect_list)]]
      return(venn_list)
    }
    venn_list = list()
    venn_list[["all"]] = lapply(names(deseq_results_sig),function(x) rownames(deseq_results_sig[[x]]))
    venn_list[["up"]] = lapply(names(deseq_results_sig),function(x) rownames(deseq_results_sig[[x]])[deseq_results_sig[[x]][["log2FoldChange"]] > 0])
    venn_list[["down"]] = lapply(names(deseq_results_sig),function(x) rownames(deseq_results_sig[[x]])[deseq_results_sig[[x]][["log2FoldChange"]] < 0])
    names(venn_list[["all"]]) = names(deseq_results_sig)
    names(venn_list[["up"]]) = names(deseq_results_sig)
    names(venn_list[["down"]]) = names(deseq_results_sig)
    lapply(names(venn_list),function(venn_name){
      lapply(names(venn_list[[venn_name]]),function(venn_write){
        write.table(x = venn_list[[venn_name]][[venn_write]],
                  file = paste0(rlog_vst,"/Venn_Significant_DEGs/Venn_Significant_DEGs_",venn_write,"_",venn_name,"_",rlog_vst,"_",genes_isoforms,".txt"),
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = FALSE)
      })
      if(length(venn_list[[venn_name]]) > 5){
      venn_sub = cut(1:length(venn_list[[venn_name]]),ceiling(length(venn_list[[venn_name]])/5), labels = FALSE)
      for(i in unique(venn_sub)){
      temp_venn = venn.diagram(venn_list[[venn_name]][venn_sub %in% i],
                               filename = NULL,
                               imagetype = "png",
                               fill = brewer.pal(brewer.pal.info["Set1","maxcolors"],"Set1")[1:sum(venn_sub %in% i)],
                               alpha = rep(0.25,sum(venn_sub %in% i)),
                               lwd = rep(0,sum(venn_sub %in% i)),
                               force.unique = TRUE,
                               ext.percent = 0.01,
                               cex = 1.5,
                               cat.cex = 1.5
      )
      png(filename = paste0("rlog/Venn_diagram_",venn_name,"_",Experiment_name,"_",i,".png"),width = 1080,height = 720,units = "px")
      grid.newpage()
      pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
      grid.draw(temp_venn)
      while(!is.null(dev.list())){dev.off()}
      }
    }else{
      temp_venn = venn.diagram(venn_list[[venn_name]],
                               filename = NULL,
                               imagetype = "png",
                               fill = brewer.pal(brewer.pal.info["Set1","maxcolors"],"Set1")[1:length(deseq_results_sig)],
                               alpha = rep(0.25,length(deseq_results_sig)),
                               lwd = rep(0,length(deseq_results_sig)),
                               force.unique = TRUE,
                               ext.percent = 0.01,
                               cex = 1.5,
                               cat.cex = 1.5
                               )
      png(filename = paste0("rlog/Venn_diagram_",venn_name,"_",Experiment_name,".png"),width = 1080,height = 720,units = "px")
      grid.newpage()
      pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
      grid.draw(temp_venn)
      while(!is.null(dev.list())){dev.off()}
    }
    })
    unlink("VennDiagram*.log")
  }
    if("Variable heatmap and report" %in% section){
  #####    Create variable distance heat map    ######
  cat("#####    Creating variable distance heat map    ######\n")
  colors=colorRampPalette(brewer.pal(9,"Blues"))(255)
  pheatmap(cor(data_set_matrix),filename = paste0(rlog_vst,"/Heatmaps/Heatmap_variables_",rlog_vst,"_",genes_isoforms,".png"),col = colors,main = "Variable heatmap",annotation_colors = color_select)

  #####    Create reports    ######
  cat("#####    Creating reports    ######\n")
  time_start=Sys.time()
  for(i in factors){
  DESeq_report = HTMLReport(shortName = paste0(rlog_vst,"_",genes_isoforms,"_",i,"_reports"),
                          title = paste0(rlog_vst,"_",genes_isoforms,"_",i,"_reports"),
                          reportDirectory = paste0(rlog_vst,"/",i,"_report"))
  publish(object = data_set_DESeq,DESeq_report,pvalueCuttoff=alpha,
          factor = colData(data_set_DESeq)[[i]],
          reportDir = paste0(rlog_vst,"/reports"),make.plots = TRUE)
  finish(DESeq_report)
  rm(DESeq_report)
  }
  
  if(length(grep(pattern = ":",design_formula))){
    for(int in grep(pattern = "_vs_",resultsNames(data_set_DESeq)[-1],value = TRUE,invert = TRUE)){
    data_set_DESeq_int = data_set_DESeq[rownames(deseq_results_sig[[int]]),]
      for(i in factors){
      DESeq_report = HTMLReport(shortName = paste0(rlog_vst,"_",genes_isoforms,"_",int,"_reports"),
                                title = paste0(rlog_vst,"_",genes_isoforms,"_",int,"_reports"),
                                reportDirectory = paste0(rlog_vst,"/",int,"_report"))
      publish(object = data_set_DESeq_int,
              publicationType = DESeq_report,
              pvalueCuttoff = alpha,
              factor = colData(data_set_DESeq)[[i]],
              reportDir = paste0(rlog_vst,"/reports"),make.plots = TRUE)
      finish(DESeq_report)
      rm(DESeq_report)
    }
    }
    }
  cat("Reports completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }
    if(exists("PCA_data")){
      save(list = c("PCA_data","experimental_design","PCA_3D","rlog_vst","genes_isoforms","wd","Experiment_name"),file = "PCA_data.RData")
    #   if(system("which xvfb-run",ignore.stdout = TRUE) == 0 & system("which ffmpeg",ignore.stdout = TRUE) == 0 & system("which 3D_PCA_run",ignore.stdout = TRUE) == 0){
    #     try(system(paste0("xvfb-run 3D_PCA_run ",paste("--factor ",factors,sep = "",collapse = " ")," --wd ",getwd())))
    #   }else{
    #     cat("Install ",paste(c("xvfb","ffmpeg","3D_PCA")[as.logical(c(system("which xvfb-run",ignore.stdout = TRUE),
    #                                                                   system("which ffmpeg", ignore.stdout = TRUE),
    #                                                                   system("which 3D_PCA",ignore.stdout = TRUE)))],
    #                          sep = "", collapse = " | "),"\n",sep = "")
    #     cat("For 3D PCA movies run: xvfb-run 3D_PCA_run --factor your_factor --wd working_dir\n",sep = "")
    #   }
    }
      save.image(paste0(Experiment_name,"_",rlog_vst,"_",genes_isoforms,"_","final_data.RData"))
  }

cat("#####    End of DEseq analysis    #####\n")

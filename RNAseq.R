#####    RNA-seq analysis (RSEM / counts)    #####

packages = c("DESeq2","ggplot2","ggrepel","gplots","RColorBrewer","BiocParallel","tximport","readr",
             "pheatmap","goseq","rstudioapi","ReportingTools","factoextra","vegan","rgl","ape","cluster","data.table",
             "parallel","doParallel","RCurl","devtools","GenomicFeatures","apeglm","R.utils","VennDiagram","wordcloud",
             "tm","topGO","Rgraphviz","NOISeq","gprofiler2","jsonlite","bcbioRNASeq")

if(!"bcbioRNASeq" %in% rownames(installed.packages())){
  BiocManager::install("MultiAssayExperiment")
  BiocManager::install("SingleCellExperiment")
  install.packages(
    pkgs = "bcbioRNASeq",
    repos = c("https://r.acidgenomics.com",getOption("repos"))
  )
}

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

if(!exists("init_params")){
  init_params = list() 
}

###### Cluster commands ######
if(!interactive()){
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  must_args = c("wd","name","process")
  if(!all(must_args %in% names(args))){
    help = matrix(data = c("--rdata    ","RData file path",
                           "--wd","Working directory path",
                           "--name","Experiment name",
                           "--input","Aligner output files path (RSEM, STAR or kallisto)",
                           "--mapper","One of RSEM, Kallisto, Salmon, HTseq-count, Counts",
                           "--design","DESeq2 design formula (e.g. ~Genotype+treatment)",
                           "--exp","Experimental design file (Control first)",
                           "--reduced","DESeq2 reduced design formula (Default = ~1)",
                           "--process","Process: ",
                           "","1 = Run settings",
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
                           "","13 = Venn diagram",
                           "","14 = topGO analysis",
                           "","15 = Word cloud",
                           "","16 = Variable heatmap and report",
                           "","(e.g. 2,5,6,7,8,9,10,11,12)",
                           "--remove_isoforms","Remove isoform suffix from gene IDs",
                           "--isoforms","Use isoforms instead of genes",
                           "--vst","Use VST instead of rlog",
                           "--lengths","File with extention .gtf, .rds, .RData, or a tab delimited file",
                           "--tx2gene","A transcript to gene mapping file (Kallisto/Salmon)",
                           "--alpha","Cut off for p values",
                           "--FDR","Cut off for FDR",
                           "--k","Number of clusters for K-means clustering",
                           "--seed","Seed value for random number generator",
                           "--heatmap_no_clust","Cluster rows in heatmaps",
                           "--GO_file","Path to GO annotation file",
                           "--noiseq","Perform NOISeq correction (ARSyNseq: counts | proportion | FDR)",
                           "--venn_go","0 = None, 1 = Up-regulated, 2 = Down-regulated, 3 = Both (Default = 0; Multiple input possible)",
                           "--ensembl","Select Ensembl species ID (enter 0 for options)",
                           "--t","Number of compute threads",
                           "--arg","Additional R arguments (multiple arguments in separate flags)")
                  ,ncol = 2,byrow = TRUE)
    prmatrix(help,quote = FALSE,rowlab = rep("",nrow(help)),collab = rep("",2))
    stop(paste0("Missing command line input --> ",paste(must_args[!must_args %in% names(args)],collapse = " | ")), call. = TRUE)
    if(sum(as.numeric(unlist(strsplit(args[["process"]],split = ","))) %in% c(2,3,4)) > 1){
      stop(paste0("Select only one DESeq2 analysis (2,3 or 4)"),call. = TRUE)
    }
  }
  if("rdata" %in% names(args)){
    cat("Loading RData file: ",args[["rdata"]],"\n",sep = "")
    load(args[["rdata"]])
    args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  }
  init_params[["wd"]] = path.expand(args[["wd"]])
  setwd(init_params[["wd"]])
  init_params[["Experiment_name"]] = args[["name"]]
  if("input" %in% names(args)){
    init_params[["input_path"]] = path.expand(args[["input"]])
  }
  
  if("exp" %in% names(args)){
    init_params[["exp_design"]] = path.expand(args[["exp"]])
  }
  
  if("mapper" %in% names(args)){
    init_params[["Mapper"]] = grep(pattern = args[["mapper"]],x = c("RSEM","Kallisto","Salmon","HTseq-count","Counts"),ignore.case = TRUE,value = TRUE)
  }
  init_params[["heatmap_row_clust"]] = !"heatmap_no_clust" %in% names(args)
  if("noiseq" %in% names(args)){
    init_params[["NOISeq_correction"]] = TRUE
  }else{
    if("NOISeq_correction" %in% names(init_params)){
      init_params[["NOISeq_correction"]] = FALSE
    }
  }
  
  if("isoforms" %in% names(args)){
    init_params[["genes_isoforms"]] = "isoforms"
  }else{
    init_params[["genes_isoforms"]] = "genes"
  }
  
  if("alpha" %in% names(args)){
    init_params[["alpha"]] = as.numeric(args[["alpha"]])
  }else{
    init_params[["alpha"]] = 0.05
  }
  
  if("FDR" %in% names(args)){
    init_params[["FDR"]] = as.numeric(args[["FDR"]])
  }else{
    init_params[["FDR"]] = 0.05
  }
  
  if("lengths" %in% names(args)){
    init_params[["lengths"]] = path.expand(args[["lengths"]])
  }
  
  if("VST" %in% names(args)){
    init_params[["rlog_vst"]] = "VST"
  }else{
    init_params[["rlog_vst"]] = "rlog"
  }
  
  if("design_formula" %in% names(args)){
    init_params[["design_formula"]] = formula(args[["design"]])
  }
  if("reduced_formula" %in% names(args)){
    if("reduced" %in% names(args)){
      init_params[["reduced_formula"]] = formula(args[["reduced"]])
    }else{
      init_params[["reduced_formula"]] = formula(~1)
    }
  }
  
  if("tx2gene" %in% names(args)){
    init_params[["tx2gene"]] = path.expand(args[["tx2gene"]])
  }
  
  init_params[["section"]] = c("Run settings",
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
              "Venn diagram",
              "topGO analysis",
              "Wordcloud",
              "Variable heatmap and report")
  
  init_params[["section"]] = init_params[["section"]][as.numeric(unlist(strsplit(args[["process"]],split = ",")))]
  if(!"Optimize K-means" %in% init_params[["section"]] & "K-means clustering" %in% init_params[["section"]]){
    if(!"k" %in% names(args)){
      stop(paste0("Error: Please specify the number of k clusters with --k OR K cluster optimization with --process 8"),call. = TRUE)
    }
    init_params[["k_clusters"]] = as.numeric(args[["k"]])
    cat("Number of K clusters: ",init_params[["k_clusters"]],"\n",sep = "")
  }else{
    cat("Performing K cluster optimization\n",sep = "")
  }
  if("GO_file" %in% names(args)){
    init_params[["GO_file"]] = path.expand(args[["GO_file"]])
  }
  
  if("ensembl" %in% names(args)){
    if(args[["ensembl"]] == 0){
      Ensembl_species = fromJSON("https://biit.cs.ut.ee/gprofiler/api/util/organisms_list/")
      prmatrix(matrix(data = Ensembl_species[order(Ensembl_species$display_name),"display_name"]),quote = FALSE,rowlab = as.character(1:nrow(Ensembl_species)),collab = rep("",2))
      stop("Select a number from the list and re-run with '--ensembl #'",call. = TRUE)
    }else{
      init_params[["Ensembl_species_ID"]] = args[["ensembl"]]
    }
  }
  
  
  if("remove_isoforms" %in% names(args)){
    init_params[["remove_isoforms"]] = TRUE
    init_params[["isoform_pattern"]] = unname(args[["remove_isoforms"]])
  }else{
    init_params[["remove_isoforms"]] = FALSE
  }
  if("seed" %in% names(args)){
    init_params[["random_seed"]] = as.numeric(args[["seed"]])
  }else{
    init_params[["random_seed"]] = 123
  }
  if("venn_go" %in% names(args)){
    init_params[["venn_GO"]] = as.numeric(unlist(strsplit(args[["venn_go"]],split = ",")))
  }else{
    if("venn_GO" %in% names(init_params)){
      init_params[["venn_GO"]] = 0
    }
  }
  if("threads" %in% names(args)){
    init_params[["threads"]] = as.numeric(args[["t"]])
  }else{
    init_params[["threads"]] = 1
  }

  ##### Additional arguments #####
  invisible(
    sapply(
      which(names(args) %in% "arg"),
      function(x){
        eval(parse(text = args[[x]]),envir = .GlobalEnv)
      }
      )
    )

  cat("Working directory: ",init_params[["wd"]],"\n",sep = "")
  cat("Experiment name: ",init_params[["Experiment_name"]],"\n",sep = "")
  cat("Number of threads: ",init_params[["threads"]],"\n",sep = "")
  cat("Running sections: ",paste(init_params[["section"]],collapse = " | "),"\n",sep = "")

  ##### Register threads #####
  if(.Platform$OS.type == "unix"){
    register(BPPARAM = MulticoreParam(workers = init_params[["threads"]]))
  }else{
    register(BPPARAM = SnowParam(workers = init_params[["threads"]]))
  }

  ##### RData output ######
  .classes = NULL
  for(.obj in ls()){
    suppressWarnings({.classes[.obj] = class(get(.obj))})
  }
  prmatrix(matrix(data = c(ls(),.classes),nrow = length(ls()),ncol = 2),quote = FALSE,rowlab = rep("",length(ls())),collab = rep("",2))
  rm(.classes,.obj)

  }else{
    
    ### Interactive

init_params[["threads"]] = detectCores()

##### Register threads #####
if(.Platform$OS.type == "unix"){
  register(BPPARAM = MulticoreParam(workers = init_params[["threads"]]))
}else{
  register(BPPARAM = SnowParam(workers = init_params[["threads"]]))
}

init_params[["section"]] = select.list(choices = c("Run settings",
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
if("Run settings" %in% init_params[["section"]]){
  ###### Set working directory ######
  cat("#####   Select working directory   #####\n")
  init_params[["wd"]] = path.expand(rstudioapi::selectDirectory(caption = "Choose working directory:"))
  setwd(init_params[["wd"]])
  
  init_params[["Experiment_name"]] = as.character(readline(prompt = "Experiment name: "))
  
  # Get mapper output and factor number
  init_params[["Mapper"]] = select.list(choices = c("RSEM","Kallisto","Salmon","HTseq-count","Counts"),multiple = FALSE,title = "Select Mapper",graphics = TRUE)
  init_params[["input_path"]] = path.expand(rstudioapi::selectDirectory(caption = "Choose mapping output directory: ",path = init_params[["wd"]]))

  if(select.list(choices = c("No","Yes"),multiple = FALSE,title = "Remove isoforms for annotation?",graphics = TRUE) == "Yes"){
    init_params[["remove_isoforms"]] = TRUE
    init_params[["isoform_pattern"]] = readline(prompt = "Input regular expression for isoform removal: ")
  }else{
    init_params[["remove_isoforms"]] = FALSE
  }
  if(select.list(choices = c("Yes","No"),multiple = FALSE,title = "Get gene lengths?",graphics = TRUE) == "Yes"){
    init_params[["lengths"]] = file.choose()
  }

  if(init_params[["Mapper"]] %in% c("Kallisto","Salmon")){
    init_params[["tx2gene"]] = selectFile(caption = "Select transcript to gene file",path = init_params[["wd"]])
  }

  if(init_params[["Mapper"]] %in% "Counts"){
    count_genes = as.numeric(readline(prompt = "Select gene column: "))
    if(!"lengths" %in% names(init_params)){
      count_length = as.numeric(readline(prompt = "Select length column: "))
    }
    count_column = as.numeric(eval(expr = parse(text = readline(prompt = "Select counts column (set range e.g. 4:7): "))))
  }

  # Create experimental design
  cat("Text file: column 1 sample names with title, column 2 file path with title, and the first row of the rest of the columns as factor names\n")
  init_params[["exp_design"]] = selectFile(caption = "Experimental design file",path = init_params[["wd"]])
  
  # Select variable reference
  cat("Independent design: ~",paste(colnames(read.table(file = init_params[["exp_design"]],header = TRUE,row.names = 1,sep = "\t"))[-1],collapse = " + "),"\n",sep = "")
  cat("Interaction term = Factor1:Factor2 \n",sep = "")
  init_params[["design_formula"]] = formula(readline(prompt = "Select main design formula: "))
  cat("Reduced design, such as batch factor. No design = ~1","\n",sep = "")
  init_params[["reduced_formula"]] = formula(readline(prompt = "Select reduced design formula: "))
  if(select.list(choices = c("Yes","No"),multiple = FALSE,title = "Cluster rows in heatmaps?",graphics = TRUE)=="Yes"){
    init_params[["heatmap_row_clust"]] = TRUE
  }else{
    init_params[["heatmap_row_clust"]] = FALSE
    }
init_params[["alpha"]] = as.numeric(readline(prompt = "Choose alpha cutoff: "))
init_params[["FDR"]] = as.numeric(readline(prompt = "Choose FDR cutoff: "))
init_params[["rlog_vst"]] = select.list(choices = c("rlog","VST"),multiple = FALSE,title = "Select transformation",graphics = TRUE)
init_params[["NOISeq_correction"]] = select.list(choices = c("No","Yes"),multiple = FALSE,title = "NOISeq correction?",graphics = TRUE) == "Yes"
init_params[["random_seed"]] = as.numeric(readline(prompt = "Select seed value for random number generator: "))
init_params[["venn_GO"]] = list("None" = 0,"Only Up-regulated" = 2,"Only Down-regulated" = 3,"Both" = 1,"All changes" = 1:3)[[select.list(choices = c("None","Only Up-regulated","Only Down-regulated","Both","All changes"),multiple = FALSE,title = "GO enrichment for Venn subsets",graphics = TRUE)]]

if(all(!"Optimize K-means" %in% init_params[["section"]],"K-means clustering" %in% init_params[["section"]])){
  init_params[["k_clusters"]] = as.numeric(readline(prompt = "Select number of clusters (k): "))
  }

GO_select = select.list(choices = c("Yes","From Ensembl","No"),multiple = FALSE,title = "Load GO annotations?",graphics = TRUE)
if(GO_select == "Yes"){
  init_params[["GO_file"]] = selectFile(caption = "Choose GO annotation file:",path = init_params[["wd"]])
}
if(GO_select == "From Ensembl"){
  Ensembl_species = fromJSON("https://biit.cs.ut.ee/gprofiler/api/util/organisms_list/")
  Ensembl_species = Ensembl_species[order(Ensembl_species$display_name),]
  init_params[["Ensembl_species_ID"]] = Ensembl_species[Ensembl_species$display_name %in% select.list(choices = Ensembl_species[,"display_name"],
                                                                multiple = FALSE),"id"]
}
    save.image(paste0(init_params[["Experiment_name"]],".RData"))
}
  }

if("Run settings" %in% init_params[["section"]]){
  
  # Experimental design
  experimental_design = read.table(file = init_params[["exp_design"]],header = TRUE,sep = "\t",row.names = 1)
  if(!all(experimental_design[[1]] %in% list.files(path = init_params[["input_path"]],full.names = TRUE))){
    cat("Changing file paths to selected input path\n")
    experimental_design[[1]] = file.path(init_params[["input_path"]],basename(experimental_design[[1]]))
    if(!all(experimental_design[[1]] %in% list.files(path = init_params[["input_path"]],full.names = TRUE))){
      stop("Error: file names in experimental design file must match file names (multiple files) / sample names (one file counts)\n",call. = TRUE)
    }}
  file_names = experimental_design[[1]]
  experimental_design = experimental_design[,-1,drop = FALSE]
  factors = colnames(experimental_design)
  for(i in colnames(experimental_design)){
    experimental_design[[i]] = relevel(as.factor(experimental_design[[i]]),ref = experimental_design[[i]][1])
  }
  write.table(x = cbind(Sample_name = rownames(experimental_design),File_name = file_names,experimental_design),file = "Experimental_design.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  
  if(init_params[["Mapper"]] %in% "RSEM"){
    init_params[["genes_isoforms"]] = if(grepl(pattern = "[.]genes[.]results",file_names[1])){
      "genes"
    }else{
      "isoforms"
    }
  }
  
  # Lengths
  if("lengths" %in% names(init_params)){
    # GTF
    if(grepl(pattern = "gtf",sub(pattern = ".+[.]",replacement = "",init_params[["lengths"]]),ignore.case = TRUE)){
    gtf = read.table(file = file.choose(),header = FALSE,sep = "\t",comment.char = "#")
    gtf = gtf[grep(pattern = "transcript|gene",gtf[,3],ignore.case = TRUE),]
    if("genes" %in% init_params[["genes_isoforms"]]){
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
    
    # RDS
    if(grepl(pattern = "rds",sub(pattern = ".+[.]",replacement = "",init_params[["lengths"]]),ignore.case = TRUE)){
      gene_length = readRDS(init_params[["lengths"]])
    }
    
    # RDATA
    if(grepl(pattern = "rdata",sub(pattern = ".+[.]",replacement = "",init_params[["lengths"]]),ignore.case = TRUE)){
      gene_length = load(init_params[["lengths"]])
    }
    
    # TAB
    if(!grepl(pattern = "gtf|rds|rdata",sub(pattern = ".+[.]",replacement = "",init_params[["lengths"]]),ignore.case = TRUE)){
      gene_length = read.table(file = init_params[["lengths"]],header = FALSE,sep = "\t",row.names = 1)
      gene_length[[1]] = as.numeric(gene_length[[1]])
      colnames(gene_length) = "length"
    }
  }
  
  # tx2gene
  if("tx2gene" %in% names(init_params)){
    tx2gene = read.table(init_params[["tx2gene"]])
  }
  
  
##### Load GO annotation #####
if("GO_file" %in% names(init_params)){
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
  geneGO = delim_fun(init_params[["GO_file"]])
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
      try(system(paste0("ffmpeg -y -hide_banner -loglevel warning -r 60 -y -i ",init_params[["wd"]],"/animation_merge/frame-%03d.png ./",init_params[["rlog_vst"]],"/PCA/PCA_3D_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",run_factor,".mp4")))
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
          try(system(paste0("ffmpeg -y -hide_banner -loglevel warning -r 60 -y -i ",init_params[["wd"]],"/animation_merge/frame-%03d.png ./",init_params[["rlog_vst"]],"/PCA/PCA_3D_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",PCA_comp,".mp4")))
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
          try(system(paste0("ffmpeg -y -hide_banner -loglevel warning -r 60 -y -i ",init_params[["wd"]],"/animation_merge/frame-%03d.png ./",init_params[["rlog_vst"]],"/PCA/PCA_3D_DEGs_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",PCA_comp,".mp4")))
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
  if("DESeq2" %in% init_params[["section"]] | "LRT" %in% init_params[["section"]]){
  cat("#####    Importing data    ######\n", sep = "")
    if("input_path" %in% names(init_params)){
      file_names = file.path(init_params[["input_path"]],basename(file_names))
    }
  if(init_params[["Mapper"]] %in% "RSEM"){
  names(file_names) = rownames(experimental_design)
  if(init_params[["genes_isoforms"]] == "genes"){
    RSEM_counts = tximport(files = file_names,type = "rsem",txIn = FALSE,txOut = FALSE)
  }
  if(init_params[["genes_isoforms"]] == "isoforms"){
    RSEM_counts = tximport(files = file_names,type = "rsem",txIn = TRUE,txOut = TRUE)
  }
  RSEM_counts$length[RSEM_counts$length == 0] = 1
  data_set = DESeqDataSetFromTximport(txi = RSEM_counts,colData = experimental_design,design = init_params[["design_formula"]])
  gene_length = RSEM_counts$length
  count_table = assay(data_set)
  }

    if(init_params[["Mapper"]] %in% "Kallisto"){
      names(file_names) = rownames(experimental_design)
      kallisto_counts = tximport(files = file_names,type = "kallisto", tx2gene = tx2gene,countsFromAbundance = "lengthScaledTPM")
      kallisto_counts$length[kallisto_counts$length == 0] = 1
      TPM = kallisto_counts$abundance
      data_set = DESeqDataSetFromTximport(txi = kallisto_counts,colData = experimental_design,design = init_params[["design_formula"]])
      gene_length = kallisto_counts$length
      count_table = assay(data_set)
    }

    if(init_params[["Mapper"]] %in% "Salmon"){
      names(file_names) = rownames(experimental_design)
      salmon_counts = tximport(files = file_names,
                           type = "salmon",
                           txOut = TRUE,
                           countsFromAbundance = "lengthScaledTPM")
      salmon_counts$length[salmon_counts$length == 0] = 1
      TPM = salmon_counts$abundance
      data_set = DESeqDataSetFromTximport(txi = salmon_counts,colData = experimental_design,design = init_params[["design_formula"]])
      gene_length = salmon_counts$length
      count_table = assay(data_set)
    }

  if(init_params[["Mapper"]] %in% "HTseq-count"){
    sampleTable = data.frame(sampleName = rownames(experimental_design),
                             filename = basename(file_names))
    sampleTable = cbind(sampleTable,experimental_design)
    data_set = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = init_params[["input_path"]],design = init_params[["design_formula"]])
    count_table = assay(data_set)
    count_table = count_table[!grepl(pattern = "^N_unmapped$|^N_multimapping$|^N_noFeature$|^N_ambiguous$",rownames(count_table)),]
  }

  if(init_params[["Mapper"]] %in% "Counts"){
    for(i in file_names){
      cat("Reading file: ",i,"\n", sep = "")
      if(i == file_names[1]){
        count_table = read.table(file = i,header = TRUE,sep = "\t")[c(count_genes,count_column)]
        rownames(count_table) = count_table[,count_genes]
        if("lengths" %in% names(init_params)){
          gene_length = read.table(file = i,header = TRUE,sep = "\t")[count_length]
        if(length(count_column) > 1){
          gene_length[,2:length(count_column)] = gene_length[[1]]
        }
        rownames(gene_length) = count_table[,count_genes]
        }
        
        count_table = count_table[,-count_genes,drop = FALSE]
      }else{
        count_table = cbind(count_table,read.table(file = i,header = TRUE,sep = "\t")[count_column])
        if("lengths" %in% names(init_params)){
          gene_length = cbind(gene_length,read.table(file = i,header = TRUE,sep = "\t")[count_length])
        }
      }
    }
        colnames(count_table) = rownames(experimental_design)

  }
    if(init_params[["NOISeq_correction"]]){
      count_table = round(ARSyNseq(data = readData(data = filtered.data(dataset = count_table,
                                                                        depth = colSums(count_table),
                                                                        p.adj = "fdr",
                                                                        method = 3,
                                                                        norm = FALSE,
                                                                        factor = factors[1]),
                                                   factors = data.frame(experimental_design,
                                                                        Replicates = apply(experimental_design,MARGIN = 1,function(x) paste(x,collapse = "."))),
                                                   length = as.matrix(gene_length)[,1]),
                                   factor = "Replicates",
                                   norm = "n",
                                   batch = FALSE,
                                   logtransf = FALSE)@assayData$exprs)
      gene_length = gene_length[rownames(count_table),,drop = FALSE]
      data_set = DESeqDataSetFromMatrix(countData = count_table,
                                        colData = experimental_design,
                                        design = init_params[["design_formula"]],
                                        tidy = FALSE,
                                        ignoreRank = FALSE)
    }

      if(!init_params[["Mapper"]] %in% c("Kallisto","Salmon")){
      TPM = if(ncol(count_table) == ncol(gene_length)){
        t(t(count_table/gene_length) * 1e6 / colSums(count_table/gene_length))
      }else{
        t(t(count_table/gene_length[,1]) * 1e6 / colSums(count_table/gene_length[,1]))
      }
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

    if("Ensembl_species_ID" %in% names(init_params)){
      geneGO = gconvert(rownames(count_table),organism = init_params[["Ensembl_species_ID"]],target = "GO")[,c("input","target")]
      gene_names = geneGO$input
      geneGO = lapply(geneGO$input,function(x) geneGO[geneGO$input %in% x,"target"])
      names(geneGO) = gene_names
    }
    
    if("species_ID" %in% names(init_params) || "GO_file" %in% names(init_params)){
      cat("##### Downloading GO terms #####\n")
      GO_terms = AnnotationDbi::select(GO.db,keys = keys(GO.db),columns = columns(GO.db))
      GO_table = data.frame(Gene = rep(names(lengths(geneGO)),lengths(geneGO)),
                            GOID = unlist(geneGO,use.names = FALSE))
      GO_table[["Term"]] = GO_terms[match(GO_table$GOID,table = GO_terms$GOID),"TERM"]
      GO_table[["Ontology"]] = GO_terms[match(GO_table$GOID,table = GO_terms$GOID),"ONTOLOGY"]
      GO_table[["Definition"]] = GO_terms[match(GO_table$GOID,table = GO_terms$GOID),"DEFINITION"]
      write.table(x = GO_table,file = paste0("GO_table_info_",init_params[["Experiment_name"]],".txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
      write.table(x = GO_table[,c("Gene","GOID")],file = paste0("GO_table_",init_params[["Experiment_name"]],".txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
      write.table(x = c("(species=Rbiotools)(type=Biological Process)(curator=GO)",paste(GO_table[GO_table$Ontology %in% "BP","Gene"],gsub(x = GO_table[GO_table$Ontology %in% "BP","GOID"],pattern = "GO:",replacement = ""),sep = " = ")),file = paste0("BinGO_BP_",init_params[["Experiment_name"]],".txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
      write.table(x = c("(species=Rbiotools)(type=Cellular Component)(curator=GO)",paste(GO_table[GO_table$Ontology %in% "CC","Gene"],gsub(x = GO_table[GO_table$Ontology %in% "CC","GOID"],pattern = "GO:",replacement = ""),sep = " = ")),file = paste0("BinGO_CC_",init_params[["Experiment_name"]],".txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
      write.table(x = c("(species=Rbiotools)(type=Molecular Function)(curator=GO)",paste(GO_table[GO_table$Ontology %in% "MF","Gene"],gsub(x = GO_table[GO_table$Ontology %in% "MF","GOID"],pattern = "GO:",replacement = ""),sep = " = ")),file = paste0("BinGO_MF_",init_params[["Experiment_name"]],".txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
    }
    
  ##### Run DESeq2 #####
  cat("#####    Running DESeq2   ######\n")
  time_start=Sys.time()
  if("DESeq2" %in% init_params[["section"]]){
    data_set_DESeq = DESeq(data_set,parallel = TRUE)
    cat("DESeq ",init_params[["genes_isoforms"]]," completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }
  if("LRT" %in% init_params[["section"]]){
    data_set_DESeq = DESeq(data_set,test = "LRT",reduced = init_params[["reduced_formula"]],parallel = TRUE)
    cat("LRT ",init_params[["genes_isoforms"]]," completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }
    save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_DESeq2_data.RData"))
  }

#####   Transform counts   #####
  if("Transformation" %in% init_params[["section"]]){
  if(init_params[["rlog_vst"]] == "rlog"){
    cat("#####    rLog transforming ",init_params[["genes_isoforms"]]," counts    ######\n", sep = "")
    time_start=Sys.time()
    data_set_transform = rlog(data_set_DESeq,blind = FALSE)
    cat("rLog transforming ",init_params[["genes_isoforms"]]," counts completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }
  if(init_params[["rlog_vst"]] == "VST"){
    cat("#####    VST transforming ",init_params[["genes_isoforms"]]," counts    ######\n", sep = "")
    time_start=Sys.time()
    data_set_transform = varianceStabilizingTransformation(data_set_DESeq,blind = FALSE)
    cat("VST transforming ",init_params[["genes_isoforms"]]," counts completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }

  # Output for WGCNA
        data_vst = if(init_params[["rlog_vst"]] == "rlog"){
          assay(varianceStabilizingTransformation(data_set_DESeq,blind = FALSE))
        }else{
          assay(data_set_transform)
        }
        data_rlog = if(init_params[["rlog_vst"]] == "rlog"){
          assay(data_set_transform)
        }else{
          assay(rlog(data_set_DESeq,blind = FALSE))
        }
        save(list = c("data_vst","data_rlog","experimental_design"), file = "WGCNA_data.RData")
        rm(data_vst,data_rlog,coldata)

  save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_","transformed_data.RData"))
  }

  if(any(c("Heatmaps","Dispersion estimates, PCAs and PCoAs","Optimize K-means","K-means clustering","Hierarchial clustering","MA-plot","DEG heatmap","goseq GO analysis","topGO analysis","Wordcloud","Venn diagram","Variable heatmap and report") %in% init_params[["section"]])){
    #####   Data matrix   #####
    dir.create(init_params[["rlog_vst"]],showWarnings = FALSE)
    dir.create(paste0(init_params[["rlog_vst"]],"/Heatmaps"),showWarnings = FALSE)
    data_set_matrix = assay(data_set_transform)
    # if(nrow(data_set_matrix) > 65536){
    #   cat("Too many genes. Removing lowest ",nrow(data_set_matrix) - 65536," expressed genes\n",sep = "")
    #   data_set_matrix = data_set_matrix[names(head(sort(abs(rowSums(data_set_matrix)),decreasing = TRUE),n = 65536)),]
    # }
    data_set_matrix_scaled = scale(data_set_matrix, center=TRUE, scale=TRUE)
    save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_","analysis_data.RData"))

    if("Heatmaps" %in% init_params[["section"]]){
    cat("#####    Creating matrices and heatmaps    ######\n")
      time_start=Sys.time()
      pheatmap(
        if(nrow(data_set_matrix) > 65536){
        data_set_matrix_scaled[names(head(sort(abs(rowSums(data_set_matrix_scaled)),decreasing = TRUE),n = 65536)),]
      }else{
        data_set_matrix_scaled
      },
             show_rownames = FALSE,
             cluster_rows = init_params[["heatmap_row_clust"]],
             annotation_col = as.data.frame(colData(data_set)),
             filename = paste0(init_params[["rlog_vst"]],"/Heatmaps/Heatmap_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],".png"))
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
               filename = paste0(init_params[["rlog_vst"]],"/Heatmaps/Heatmap_gene_subset_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],".png")
               )
    }
    save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_","analysis_data.RData"))
    cat("Matrices and heatmaps completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
    }

  if("Dispersion estimates, PCAs and PCoAs" %in% init_params[["section"]]){
  #####   Dispersion estimates   #####
  time_start=Sys.time()
  cat("#####    Creating dispersion estimates    ######\n")
  png(filename = paste0("Plot_dispertion_estimates_",init_params[["genes_isoforms"]],".png"),width = 1920,height = 1080,units = "px")
  plotDispEsts(data_set_DESeq)
  while (!is.null(dev.list())){dev.off()}
  cat("Dispersion estimates completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  #####   Single factor PCA   #####
  time_start=Sys.time()
  dir.create(paste0(init_params[["rlog_vst"]],"/PCA"),showWarnings = FALSE)
  PCA_data = list()
  cat("#####    Creating PCA and factor heatmap plots    ######\n")
  for(i in factors){
      PCA_data[["Single_factor"]][[i]] = plotPCA_PC123(object = data_set_transform,intgroup=i,returnData = TRUE)
      if(!is.null(PCA_data[["Single_factor"]][[i]])){
      png(filename = paste0(init_params[["rlog_vst"]],"/PCA/PCA_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",i,".png"),width = 1920,height = 1080,units = "px")
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


      png(filename = paste0(init_params[["rlog_vst"]],"/PCA/PCA_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",i,"_PC2.png"),width = 1920,height = 1080,units = "px")
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
    png(filename = paste0(init_params[["rlog_vst"]],"/PCA/PCA_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",i,"_vs_",j,".png"),width = 1920,height = 1080,units = "px")
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

    png(filename = paste0(init_params[["rlog_vst"]],"/PCA/PCA_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",i,"_vs_",j,"_PC2.png"),width = 1920,height = 1080,units = "px")
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
    png(filename = paste0(init_params[["rlog_vst"]],"/PCA/PCoA_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",i,"_",j,".png"),width = 1920,height = 1080,units = "px")
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
  save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_","analysis_data.RData"))
  cat("PCAs and PCoAs completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }

  #######      Calculate optimal K      #######
  if("Optimize K-means" %in% init_params[["section"]]){
    time_start=Sys.time()
    cat("#####    Calculating optimal K   ######\n")
    cpucores = init_params[["threads"]]
    options("mc.cores" = cpucores)
    registerDoParallel(cpucores)
    k_clusGap = clusGapKB(x = data_set_matrix_scaled,FUNcluster = kmeans,iter.max = 50,K.max = 50,B = if(exists("k_boot_max")){k_boot_max}else{100})
    init_params[["k_clusters"]] = with(k_clusGap, maxSE(Tab[,"gap"], Tab[,"SE.sim"], method="firstSEmax"))
    png(filename = paste0("./clusGap_k_",init_params[["k_clusters"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],".png"),width = 1920,height = 1080,units = "px")
    plot(k_clusGap, frame = FALSE, xlab = "k clusters", ylab = "Gap statistic",main = paste0("k clusters vs. gap (k = ",init_params[["k_clusters"]],")"))
    abline(v = init_params[["k_clusters"]])
    while (!is.null(dev.list())){dev.off()}

    cat("Using K = ",init_params[["k_clusters"]],"\n", sep = "")
    cat("Completed K optimization in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
    save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_","analysis_data.RData"))
  }

  #######      K-means clustering      #########
  if("K-means clustering" %in% init_params[["section"]]){
    cat("#####    K-means clustering and creating matplots    ######\n")
    time_start=Sys.time()
      set.seed(init_params[["random_seed"]])
      km.res = kmeans(x = data_set_matrix_scaled,centers = init_params[["k_clusters"]], iter.max = 100)
      samples_df = setDT(as.data.frame(data_set_matrix_scaled), keep.rownames = TRUE)[]
      km.table = setDT(as.data.frame(km.res$cluster), keep.rownames = TRUE)[]
      colnames(km.table)[2]="Cluster"
      write.table(km.table,file = "kmeans_cluster_table.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
      for(i in 1:init_params[["k_clusters"]]){
        cluster_exp = samples_df[samples_df$rn %in% km.table[km.table$Cluster==i,]$rn,]
        cluster_exp = cluster_exp[,-1]
        png(filename = paste0(init_params[["rlog_vst"]],"/Matplot_Cluster_",i,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],".png"),width = 1920,height = 1080,units = "px")
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
                    filename = paste0(init_params[["rlog_vst"]],"/Heatmaps/Heatmap_kclusters","_",init_params[["k_clusters"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],".png"),
                    kmeans_k = init_params[["k_clusters"]],
                    seed_input = init_params[["random_seed"]])
      save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_","analysis_data.RData"))
      cat("Completed K-means clustering in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }

  ##### Create and plot Hierarchical clusters #####
  if("Hierarchial clustering" %in% init_params[["section"]]){
    cat("#####    Hierarchial clustering    ######\n")
    time_start=Sys.time()
    hclusters = hclust(dist(
      if(nrow(data_set_matrix) > 65536){
        data_set_matrix_scaled[names(head(sort(abs(rowSums(data_set_matrix_scaled)),decreasing = TRUE),n = 65536)),]
      }else{
        data_set_matrix
      }, method = "euclidean"), method = "complete")
    png(filename = paste0("./Hierarchical_Clustering","_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],".png"),width = 1920,height = 1080,units = "px")
    plot(hclusters, labels = FALSE, hang = -1,ann = FALSE)
    title(xlab = "Genes",ylab = "Euclidean distance")
    if("k_clusters" %in% names(init_params)){
      rh = rect.hclust(hclusters, k = init_params[["k_clusters"]], border = 1:init_params[["k_clusters"]],cluster = cutree(hclusters, k = init_params[["k_clusters"]]))
      rh_lengths = head(cumsum(c(1, lengths(rh))), -1)
      rh_y = c(rep(-1,times = init_params[["k_clusters"]]))
      for(i in 2:init_params[["k_clusters"]]){
        if(rh_lengths[[i-1]] > rh_lengths[[i]] - 50){
          rh_y[[i]] = rh_y[[i-1]] + 0.5
          }
      }
      text(x = rh_lengths+35, y = rh_y, col="black", labels=1:init_params[["k_clusters"]], font=1)
    }
    while (!is.null(dev.list())){dev.off()}
    cat("Completed Hierarchical clustering in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
    save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_","analysis_data.RData"))
  }

  #####   MA summary, MA plots, Gene heat map and GO annotation #####
  if("MA-plot" %in% init_params[["section"]]){
    time_start=Sys.time()
    for(i in resultsNames(data_set_DESeq)[-1]){
    cat("#####    Shrinking results and calculating ,MA summary and plots    ######\n")
    deseq_results_LFC = lfcShrink(data_set_DESeq,coef = i,type = "apeglm",parallel = TRUE)
    png(filename = paste0("./MA-plot_",i,"_LFC_deseq_results_",init_params[["genes_isoforms"]],"_padj.png"),width = 1920,height = 1080,units = "px")
    plotMA(deseq_results_LFC, alpha = init_params[["alpha"]], main = paste0(i), ylim = c(min(deseq_results_LFC$log2FoldChange),max(deseq_results_LFC$log2FoldChange)))
    while (!is.null(dev.list())){dev.off()}
    }
    save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_","analysis_data.RData"))
    cat("MA summary and MA plots completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }

  # if(init_params[["remove_isoforms"]]){
  #   gene_length = matrix(data = unlist(bplapply(X = unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = rownames(gene_length))),FUN = function(x, y = gene_length){
  #     return(colMaxs(gene_length[grepl(pattern = x,x = rownames(gene_length)),,drop = FALSE]))
  #   })),ncol = ncol(gene_length),dimnames = list(unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = rownames(gene_length))),colnames(gene_length)))
  #   save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_","analysis_data.RData"))
  # }

  ##### DEGs #####
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
      sink(file = paste0(init_params[["genes_isoforms"]],"_summary.txt"),append = TRUE)
      cat(i," vs. ",control_var, sep = "")
      summary(object = temp_results)
      sink()
      deseq_results[[compare_var]] = as.data.frame(temp_results)
      attr(deseq_results[[compare_var]],which = "factor") = f
      deseq_results[[compare_var]] = deseq_results[[compare_var]][order(deseq_results[[compare_var]][["padj"]]),]
      cat("Note: Removing ",sum(is.na(deseq_results[[compare_var]]))," empty genes from analysis\n", sep = "")
      deseq_results[[compare_var]] = deseq_results[[compare_var]][!is.na(deseq_results[[compare_var]][["padj"]]),]
      deseq_results_sig[[compare_var]] = deseq_results[[compare_var]][deseq_results[[compare_var]][["padj"]] < init_params[["alpha"]] & abs(deseq_results[[compare_var]][["log2FoldChange"]]) >= 1,]
      write.table(deseq_results_sig[[compare_var]], file = paste0(compare_var,"_deseq_results_",init_params[["genes_isoforms"]],"_padj.txt"),quote = FALSE,sep = "\t")
    }
    rm(temp_results)

  }

  ## Interactions
  if(length(grep(pattern = ":",init_params[["design_formula"]]))){
  for(compare_var in grep(pattern = "_vs_",resultsNames(data_set_DESeq)[-1],value = TRUE,invert = TRUE)){
    cat("#####    Getting results for interaction: ",compare_var,"   ######\n",sep = "")
    temp_results = results(data_set_DESeq, name = compare_var,parallel = TRUE)
    sink(file = paste0(init_params[["genes_isoforms"]],"_summary.txt"),append = TRUE)
    cat(compare_var,"_interaction", sep = "")
    summary(object = temp_results,alpha = init_params[["alpha"]])
    sink()
    deseq_results[[compare_var]] = as.data.frame(temp_results)
    attr(deseq_results[[compare_var]],which = "factor") = compare_var
    rm(temp_results)
    deseq_results[[compare_var]] = deseq_results[[compare_var]][order(deseq_results[[compare_var]][["padj"]]),]
    cat("Note: Removing ",sum(is.na(deseq_results[[compare_var]]))," empty genes from analysis\n", sep = "")
    deseq_results[[compare_var]] = deseq_results[[compare_var]][!is.na(deseq_results[[compare_var]][["padj"]]),]
    deseq_results_sig[[compare_var]] = deseq_results[[compare_var]][deseq_results[[compare_var]][["padj"]] < init_params[["alpha"]] & abs(deseq_results[[compare_var]][["log2FoldChange"]]) >= 1,]
    write.table(deseq_results_sig[[compare_var]], file = paste0(compare_var,"_deseq_results_",init_params[["genes_isoforms"]],"_padj.txt"),quote = FALSE,sep = "\t")
  }
  }
  save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_","analysis_data.RData"))

  dir.create(paste0(init_params[["rlog_vst"]],"/Significant_DEGs"),showWarnings = FALSE)
  dir.create(paste0(init_params[["rlog_vst"]],"/PCA"),showWarnings = FALSE)
  if(!exists("PCA_data")){
    PCA_data = list()
  }
  for(compare_var in grep(pattern = "_vs_",resultsNames(data_set_DESeq)[-1],value = TRUE)){

      #####    Create count matrix heat map    ######
      deseq_sig = rownames(deseq_results_sig[[compare_var]])
      DE_genes_sig = deseq_sig
      DE_genes_sig_up = rownames(deseq_results_sig[[compare_var]][deseq_results_sig[[compare_var]][["log2FoldChange"]] > 0,,drop = FALSE])
      DE_genes_sig_down = rownames(deseq_results_sig[[compare_var]][deseq_results_sig[[compare_var]][["log2FoldChange"]] < 0,,drop = FALSE])
      if(init_params[["remove_isoforms"]]){
        DE_genes_sig = unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = DE_genes_sig))
        DE_genes_sig_up = unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = DE_genes_sig_up))
        DE_genes_sig_down = unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = DE_genes_sig_down))
      }
      write.table(x = DE_genes_sig,
                  file = paste0(init_params[["rlog_vst"]],"/Significant_DEGs/Significant_DEGs_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_all.txt"),
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = FALSE)
      write.table(x = DE_genes_sig_up,
                  file = paste0(init_params[["rlog_vst"]],"/Significant_DEGs/Significant_DEGs_up_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_all.txt"),
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = FALSE)
      write.table(x = DE_genes_sig_down,
                  file = paste0(init_params[["rlog_vst"]],"/Significant_DEGs/Significant_DEGs_down_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_all.txt"),
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = FALSE)

      if("k_clusters" %in% names(init_params)){
        for(k in 1:k_clusters){
        DE_genes_sig = deseq_sig[deseq_sig %in% km.table[km.table$Cluster == k,]$rn]
        if(init_params[["remove_isoforms"]]){
          DE_genes_sig = unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = DE_genes_sig))
        }
        write.table(x = DE_genes_sig,
                    file = paste0(init_params[["rlog_vst"]],"/Significant_DEGs/Significant_DEGs_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,".txt"),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE)
        }
      }
      if("DEG heatmap" %in% init_params[["section"]] & length(deseq_sig) > 0){
      cat("#####    Creating DEG count matrix heat map    ######\n")
      pheatmap(data_set_matrix_scaled[deseq_sig,,drop = FALSE],
               show_rownames = FALSE,
               cluster_rows = init_params[["heatmap_row_clust"]],
               annotation_col = as.data.frame(colData(data_set)),
               filename = paste0(init_params[["rlog_vst"]],"/Heatmaps/Heatmap_significant_DEGs_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],".png")
               )
      if("k_clusters" %in% names(init_params)){
        if(length(deseq_sig) >= init_params[["k_clusters"]]){
        pheatmap_seed(data_set_matrix_scaled[deseq_sig,,drop = FALSE],
                      show_rownames = TRUE,
                      annotation_col = as.data.frame(colData(data_set)),
                      filename = paste0(init_params[["rlog_vst"]],"/Heatmaps/Heatmap_significant_DEGs_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_k",init_params[["k_clusters"]],".png"),
                      kmeans_k = init_params[["k_clusters"]],
                      seed_input = init_params[["random_seed"]])
        }}
        save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_","analysis_data.RData"))
        cat(" DEG heatmap completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
      }

      ######       PCA DEGs       ######
      if("Dispersion estimates, PCAs and PCoAs" %in% init_params[["section"]]){

        # Create Volcano plot
        volcano_data = deseq_results[[compare_var]][,c("log2FoldChange","padj")]
      volcano_data[["Expression"]] = apply(volcano_data,MARGIN = 1,function(x){
        if(abs(as.numeric(x[1])) >= 1 & as.numeric(x[2]) < 0.05){
          if(as.numeric(x[1]) > 0){
            "Up-regulated"
          }else{
            "Down-regulated"
          }
          }else{
          "Unchanged"
          }
        })

      volcano_plot = ggplot(volcano_data, aes(log2FoldChange, -log(padj,10))) +
        geom_point(aes(color = Expression,alpha = Expression), size = 3) +
        xlab(expression("log"[2]*"FC")) +
        ylab(expression("-log"[10]*"FDR")) +
        scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
        scale_alpha_manual(values = c(1,0.25,1)) +
        guides(colour = guide_legend(override.aes = list(size = 3))) +
        geom_hline(yintercept = -log(0.05,10), linetype="dashed") +
        geom_vline(xintercept = c(-1,1), linetype="dashed") +
        theme(text = element_text(size = 30))

      png(filename = paste0("Volcano_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",compare_var,".png"),width = 1920,height = 1080,units = "px")
      print(volcano_plot)
      while (!is.null(dev.list())){dev.off()}

      if(!is.null(attr(deseq_results[[compare_var]],which = "factor")) & length(deseq_sig) > 0){
      PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]] = plotPCA_PC123(object = data_set_transform[rownames(deseq_results_sig[[compare_var]]),],intgroup=attr(deseq_results[[compare_var]],which = "factor"),returnData = TRUE)
      if(!is.null(PCA_data[["DEGs"]][[attr(deseq_results[[compare_var]],which = "factor")]][[compare_var]])){
      png(filename = paste0(init_params[["rlog_vst"]],"/PCA/PCA_DEGs_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",compare_var,".png"),width = 1920,height = 1080,units = "px")
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

      png(filename = paste0(init_params[["rlog_vst"]],"/PCA/PCA_DEGs_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",compare_var,"_PC2.png"),width = 1920,height = 1080,units = "px")
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
      if("goseq GO analysis" %in% init_params[["section"]]){
        time_start=Sys.time()
        cat("#####    Starting goseq GO annotation    ######\n")
        dir.create(paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Annotation_table_sig_info"),showWarnings = FALSE, recursive = TRUE)
        dir.create(paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Annotation_table_sig"),showWarnings = FALSE)
        dir.create(paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Bootstrap"),showWarnings = FALSE)
        dir.create(paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Enrichment"),showWarnings = FALSE)
        dir.create(paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/PWF"),showWarnings = FALSE)
        dir.create(paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Wordcloud"),showWarnings = FALSE)
        for(k in if("k_clusters" %in% names(init_params)){
                  c("all",1:init_params[["k_clusters"]])
                  }else{
                  "all"}){
          if(k == "all"){
            DE_genes_sig = deseq_sig
            if(init_params[["remove_isoforms"]]){
              DE_genes_sig = unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = DE_genes_sig))
            }
            write.table(x = GO_table[GO_table[["Gene"]] %in% DE_genes_sig,],
                        file = paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Annotation_table_sig_info/Annotation_table_sig_info_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_all.txt"),
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        sep = "\t")
            write.table(x = GO_table[GO_table[["Gene"]] %in% DE_genes_sig,c("Gene","GOID")],
                        file = paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Annotation_table_sig/Annotation_table_sig_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_all.txt"),
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        sep = "\t")
          }else{
            DE_genes_sig = deseq_sig[deseq_sig %in% km.table[km.table$Cluster == k,]$rn]
            if(init_params[["remove_isoforms"]]){
              DE_genes_sig = unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = DE_genes_sig))
            }
            write.table(x = GO_table[GO_table[["Gene"]] %in% DE_genes_sig,],
                        file = paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Annotation_table_sig_info/Annotation_table_sig_info_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,".txt"),
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        sep = "\t")
            write.table(x = GO_table[GO_table[["Gene"]] %in% DE_genes_sig,c("Gene","GOID")],
                        file = paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Annotation_table_sig/Annotation_table_sig_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,".txt"),
                        quote = FALSE,
                        row.names = FALSE,
                        col.names = FALSE,
                        sep = "\t")
          }

        all_genes = rownames(deseq_results[[compare_var]])
        if(init_params[["remove_isoforms"]]){
          rownames(gene_length) = sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = rownames(gene_length))
          all_genes = unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = all_genes))
        }
        gene_vector = as.integer(all_genes %in% DE_genes_sig)

        # Create vector of gene lengths
        names(gene_vector) = all_genes
        gene_length_vector = rowMeans(gene_length[match(all_genes,rownames(gene_length)),,drop = FALSE])

        # Create Probability Weighting Function (PWF)
        png(filename = paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/PWF/PWF_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,".png"),width = 1920,height = 1080,units = "px")
        pwf = nullp(DEgenes = gene_vector,bias.data = gene_length_vector)
        while (!is.null(dev.list())){dev.off()}
        rownames(pwf) = names(gene_length_vector)
        write.table(x = pwf,file = paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/PWF/PWF_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,".txt"),quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)

        # Estimate p-values for length-corrected DE genes
        pval = goseq(pwf = pwf,gene2cat = GO_table,method = "Wallenius")
        write.table(x = pval,file = paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Enrichment/GO_enriched_table_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,".txt"),quote = FALSE,row.names = FALSE)

        # Bootstrap for comparison to p-value estimation
        GO_sample = goseq(pwf = pwf,gene2cat = GO_table,method = "Sampling",repcnt = 2000)
        png(filename = paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Bootstrap/GO_bootstrap_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,".png"),width = 1920,height = 1080,units = "px")
        plot(log10(pval[,2]), log10(GO_sample[match(pval[,1],GO_sample[,1]),2]),xlab=expression('Log'[10]*' (Wallenius p-values)'),ylab=expression('Log'[10]*' (Sampling p-values)'),xlim=c(-3,0))
        abline(0,1,col = 3,lty = 2)
        while (!is.null(dev.list())){dev.off()}

        # Export table of FDR-corrected GO annotations
        write.table(x = pval$category[p.adjust(pval$over_represented_pvalue,method="BH") < init_params[["FDR"]]],
                    file = paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Enrichment/GO_over_enriched_category_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,".txt"),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE)
        write.table(x = pval$term[p.adjust(pval$over_represented_pvalue,method="BH") < init_params[["FDR"]]],
                    file = paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Enrichment/GO_over_enriched_term_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,".txt"),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE)
        write.table(x = pval$category[p.adjust(pval$under_represented_pvalue,method="BH") < init_params[["FDR"]]],
                    file = paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Enrichment/GO_under_enriched_category_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,".txt"),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE)
        write.table(x = pval$term[p.adjust(pval$under_represented_pvalue,method="BH") < init_params[["FDR"]]],
                    file = paste0(init_params[["rlog_vst"]],"/goseq_GO_annotation/Enrichment/GO_under_enriched_term_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,".txt"),
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE)

        }
        cat("goseq GO analysis completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
      }

      ##### Venn Diagrams #####
      if("Venn diagram" %in% init_params[["section"]]){
        # Intersect function
         venn_calculate = function(deg_list){
          intersect_list = list()
          res_list = list()
          intersect_list[[1]] = deg_list
          for(i in 1:length(intersect_list[[1]])){
            attr(x = intersect_list[[1]][[i]],which = "Sample") = names(intersect_list[[1]][i])
          }
          if(length(intersect_list[[1]]) > 1){
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
          for(intersect_order in 1:(length(intersect_list)-1)){
            res_list[[intersect_order]] = intersect_list[[intersect_order]]
            for(sample_query in 1:length(intersect_list[[intersect_order]])){
              for(sample_subject in 1:length(intersect_list[[intersect_order+1]])){
                if(any(attr(intersect_list[[intersect_order]][[sample_query]],which = "Sample") %in% attr(intersect_list[[intersect_order+1]][[sample_subject]],which = "Sample"))){
                  res_list[[intersect_order]][[sample_query]] = res_list[[intersect_order]][[sample_query]][!res_list[[intersect_order]][[sample_query]] %in% intersect_list[[intersect_order+1]][[sample_subject]]]
                  attr(res_list[[intersect_order]][[sample_query]], which = "Sample") = attr(intersect_list[[intersect_order]][[sample_query]],which = "Sample")
                }
              }
            }
          }
          }
          res_list[[length(intersect_list)]] = intersect_list[[length(intersect_list)]]
          return(res_list)
        }

        venn_list = list()
        venn_list[["all"]] = lapply(names(deseq_results_sig),function(x) rownames(deseq_results_sig[[x]]))
        venn_list[["up"]] = lapply(names(deseq_results_sig),function(x) rownames(deseq_results_sig[[x]])[deseq_results_sig[[x]][["log2FoldChange"]] > 0])
        venn_list[["down"]] = lapply(names(deseq_results_sig),function(x) rownames(deseq_results_sig[[x]])[deseq_results_sig[[x]][["log2FoldChange"]] < 0])
        names(venn_list[["all"]]) = names(deseq_results_sig)
        names(venn_list[["up"]]) = names(deseq_results_sig)
        names(venn_list[["down"]]) = names(deseq_results_sig)

        venn_calc = lapply(venn_list,venn_calculate)
        names(venn_calc) = names(venn_list)


        lapply(names(venn_list),function(venn_name){

          # lapply(names(venn_list[[venn_name]]),function(venn_write){
          #   write.table(x = venn_list[[venn_name]][[venn_write]],
          #               file = paste0(init_params[["rlog_vst"]],"/Significant_DEGs/Venn_Significant_DEGs_",venn_write,"_",venn_name,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],".txt"),
          #               quote = FALSE,
          #               row.names = FALSE,
          #               col.names = FALSE)
          # })

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
              png(filename = paste0("rlog/Venn_diagram_",venn_name,"_",init_params[["Experiment_name"]],"_",i,".png"),width = 1080,height = 720,units = "px")
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
            png(filename = paste0("rlog/Venn_diagram_",venn_name,"_",init_params[["Experiment_name"]],".png"),width = 1080,height = 720,units = "px")
            grid.newpage()
            pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
            grid.draw(temp_venn)
            while(!is.null(dev.list())){dev.off()}
          }
        })
        unlink("VennDiagram*.log")
      }

      #######     topGO annotation      #########
      if("topGO analysis" %in% init_params[["section"]]){
        time_start=Sys.time()
        topGO_fun = function(DE_genes_sig,compare_var,k){
          geneList = factor(as.integer(names(geneGO) %in% DE_genes_sig),levels = c(0,1))
          names(geneList) = names(geneGO)
          for(ontology in c("BP","MF","CC")){
            sampleGOdata = new("topGOdata",
                               description = "Simple session",
                               ontology = ontology,
                               allGenes = geneList,
                               annot = annFUN.gene2GO,
                               gene2GO = geneGO)
            enrich_result = runTest(sampleGOdata, statistic = "fisher",algorithm = "classic")

            if(sum(enrich_result@score <= init_params[["FDR"]]) > 0){
              printGraph(object = sampleGOdata,
                         result = enrich_result,
                         firstSigNodes = 10,
                         fn.prefix = paste0(init_params[["rlog_vst"]],"/top_GO_annotation/topGO_graphs/",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_Best10_kcluster_",k,"_",ontology),
                         useInfo = "all",
                         pdfSW = TRUE)
              GO_DEGs_df = GenTable(object = sampleGOdata,classicFisher = enrich_result,topNodes = sum(enrich_result@score <= init_params[["FDR"]]))
              GO_DEGs_df[["Term"]] = GO_terms$TERM[match(GO_DEGs_df$GO.ID,table = GO_terms$GOID)]
              GO_DEGs_df[["Definition"]] = GO_terms$DEFINITION[match(GO_DEGs_df$GO.ID,table = GO_terms$GOID)]
              write.table(GO_DEGs_df,file = paste0(init_params[["rlog_vst"]],"/top_GO_annotation/topGO_DEGs/",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,"_",ontology,"_sig.txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

              if("Wordcloud" %in% init_params[["section"]]){
                dir.create(paste0(init_params[["rlog_vst"]],"/Wordcloud"),showWarnings = FALSE)
                cat("Creating enriched DEG Wordclouds...\n",sep = "")
                DEG_table = table(GO_DEGs_df[,"Term"])
                DEG_table = DEG_table[grep(pattern = "biological_process|cellular_component|molecular_function",x = names(DEG_table),invert = TRUE)]
                cat(format(round((3*match(k,if("k_clusters" %in% names(init_params)){c("all",1:init_params[["k_clusters"]])}else{"all"}))/(3*if("k_clusters" %in% names(init_params)){init_params[["k_clusters"]] + 1}else{1}),digits = 2)*100,nsmall = 0),"% --> Comparison: ",compare_var," | k: ",k," | Ontology: ",ontology,"           \r",sep = "")
                if(length(DEG_table) > 0){
                  png(filename = paste0(init_params[["rlog_vst"]],"/Wordcloud/Enriched_Wordcloud_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,"_",ontology,".png"),width = 1080,height = 1080,units = "px")
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
        dir.create(paste0(init_params[["rlog_vst"]],"/top_GO_annotation/topGO_graphs"),showWarnings = FALSE, recursive = TRUE)
        dir.create(paste0(init_params[["rlog_vst"]],"/top_GO_annotation/topGO_DEGs"),showWarnings = FALSE)
        if(exists("venn_calc") & any(init_params[["venn_GO"]] > 1)){
          lapply(names(venn_calc)[init_params[["venn_GO"]]],function(x){
                   lapply(venn_calc[[x]],function(y){
                     sapply(y,function(z){
                       write.table(z,file = paste0(init_params[["rlog_vst"]],"/Significant_DEGs/Venn_",x,"_",paste(attr(z,which = "Sample"),collapse = "--")),quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
                       if(init_params[["remove_isoforms"]]){
                         z = unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = z))
                       }
                       topGO_fun(DE_genes_sig = z,
                                 compare_var = paste0("Venn_",x,"_",paste(attr(z,which = "Sample"),collapse = "--")),
                                 k = "venn")
                                 })
                     })
            })
        }
        cat("#####    Starting topGO annotation    ######\n")
        for(k in if("k_clusters" %in% names(init_params)){
          c("all",1:init_params[["k_clusters"]])
        }else{
          "all"}){
          if(k == "all"){
            DE_genes_sig = deseq_sig
            if(init_params[["remove_isoforms"]]){
              DE_genes_sig = unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = DE_genes_sig))
            }
          }else{
            DE_genes_sig = deseq_sig[deseq_sig %in% km.table[km.table$Cluster == k,]$rn]
            if(init_params[["remove_isoforms"]]){
              DE_genes_sig = unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = DE_genes_sig))
            }
          }
          if(sum(names(geneGO) %in% DE_genes_sig) > 0){
            topGO_fun(DE_genes_sig = DE_genes_sig,compare_var = compare_var,k = k)
          }
        }
        # try(system(command = paste0("for i in $(find . -type f -name "*.ps"); do convert -density 800 -rotate 90 $i ${i/.ps/.png} && rm -f $i; done")))
        cat("topGO analysis completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
      }

      ##### Create Wordcloud #####
      if("Wordcloud" %in% init_params[["section"]]){
          dir.create(paste0(init_params[["rlog_vst"]],"/Wordcloud"),showWarnings = FALSE)
        cat("Creating DEG Wordclouds...\n",sep = "")
        for(k in if("k_clusters" %in% names(init_params)){
          c("all",1:init_params[["k_clusters"]])
        }else{
          "all"}){
          if(k == "all"){
            DE_genes_sig = deseq_sig
            if(init_params[["remove_isoforms"]]){
              DE_genes_sig = unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = DE_genes_sig))
            }
          }else{
            DE_genes_sig = deseq_sig[deseq_sig %in% km.table[km.table$Cluster == k,]$rn]
            if(init_params[["remove_isoforms"]]){
              DE_genes_sig = unique(sub(pattern = init_params[["isoform_pattern"]],replacement = "",x = DE_genes_sig))
            }
          }
        for(ontology in c("BP","MF","CC")){
          DEG_table = GO_table[GO_table[["Gene"]] %in% DE_genes_sig,]
          DEG_table = DEG_table[!is.na(DEG_table[["Term"]]),,drop = FALSE]
          DEG_table = table(DEG_table[DEG_table[["Ontology"]] == ontology,"Term"])
          DEG_table = DEG_table[grep(pattern = "biological_process|cellular_component|molecular_function",x = names(DEG_table),invert = TRUE)]
          cat(format(round((3*match(k,if("k_clusters" %in% names(init_params)){c("all",1:init_params[["k_clusters"]])}else{"all"}))/(3*if("k_clusters" %in% names(init_params)){init_params[["k_clusters"]] + 1}else{1}),digits = 2)*100,nsmall = 0),"% --> Comparison: ",compare_var," | k: ",k," | Ontology: ",ontology,"           \r",sep = "")
          if(length(DEG_table) > 0){
            png(filename = paste0(init_params[["rlog_vst"]],"/Wordcloud/Wordcloud_",compare_var,"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_kcluster_",k,"_",ontology,".png"),width = 1080,height = 1080,units = "px")
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
  save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_","analysis_data.RData"))
  }

    if("Variable heatmap and report" %in% init_params[["section"]]){
  #####    Create variable distance heat map    ######
  cat("#####    Creating variable distance heat map    ######\n")
  colors=colorRampPalette(brewer.pal(9,"Blues"))(255)
  pheatmap(cor(data_set_matrix),filename = paste0(init_params[["rlog_vst"]],"/Heatmaps/Heatmap_variables_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],".png"),col = colors,main = "Variable heatmap")

  #####    Create reports    ######
  cat("#####    Creating reports    ######\n")
  time_start=Sys.time()
  for(i in factors){
  DESeq_report = HTMLReport(shortName = paste0(init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",i,"_reports"),
                          title = paste0(init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",i,"_reports"),
                          reportDirectory = paste0(init_params[["rlog_vst"]],"/",i,"_report"))
  publish(object = data_set_DESeq,
          DESeq_report,
          pvalueCutoff=init_params[["alpha"]],
          factor = colData(data_set_DESeq)[[i]],
          reportDir = paste0(init_params[["rlog_vst"]],"/reports"),
          make.plots = TRUE)
  finish(DESeq_report)
  rm(DESeq_report)
  }

  if(length(grep(pattern = ":",init_params[["design_formula"]]))){
    for(int in grep(pattern = "_vs_",resultsNames(data_set_DESeq)[-1],value = TRUE,invert = TRUE)){
    data_set_DESeq_int = data_set_DESeq[rownames(deseq_results_sig[[int]]),]
      for(i in factors){
      DESeq_report = HTMLReport(shortName = paste0(init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",int,"_",i,"_reports"),
                                title = paste0(init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_",int,"_",i,"_reports"),
                                reportDirectory = paste0(init_params[["rlog_vst"]],"/",int,"_",i,"_report"))
      publish(object = data_set_DESeq_int,
              publicationType = DESeq_report,
              pvalueCutoff = init_params[["alpha"]],
              factor = colData(data_set_DESeq)[[i]],
              reportDir = paste0(init_params[["rlog_vst"]],"/reports"),
              make.plots = TRUE)
      finish(DESeq_report)
      rm(DESeq_report)
    }
    }
    }
  cat("Reports completed in ",format(round(Sys.time()-time_start,2),nsmall=2),"\n", sep = "")
  }
    if(exists("PCA_data")){
      save(list = c("PCA_data","experimental_design","PCA_3D","init_params"),file = "PCA_data.RData")
    #   if(system("which xvfb-run",ignore.stdout = TRUE) == 0 & system("which ffmpeg",ignore.stdout = TRUE) == 0 & system("which 3D_PCA_run",ignore.stdout = TRUE) == 0){
    #     try(system(paste0("xvfb-run 3D_PCA_run ",paste("--factor ",factors,sep = "",collapse = " ")," --wd ",init_params[["wd"]])))
    #   }else{
    #     cat("Install ",paste(c("xvfb","ffmpeg","3D_PCA")[as.logical(c(system("which xvfb-run",ignore.stdout = TRUE),
    #                                                                   system("which ffmpeg", ignore.stdout = TRUE),
    #                                                                   system("which 3D_PCA",ignore.stdout = TRUE)))],
    #                          sep = "", collapse = " | "),"\n",sep = "")
    #     cat("For 3D PCA movies run: xvfb-run 3D_PCA_run --factor your_factor --wd working_dir\n",sep = "")
    #   }
    }
      save.image(paste0(init_params[["Experiment_name"]],"_",init_params[["rlog_vst"]],"_",init_params[["genes_isoforms"]],"_","final_data.RData"))
  }

cat("#####    End of DEseq analysis    #####\n")

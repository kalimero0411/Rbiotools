##### ChIP analysis #####

packages=c("BiocParallel","parallel","ShortRead","RSQLite","QuasR","BSgenome","rstudioapi","ggupset",
           "Rsamtools","rtracklayer","GenomicFeatures","txdbmaker","Hmisc","Gviz","XML","mosaics",
           "ChIPseeker","clusterProfiler","ReactomePA","dada2","RMariaDB","DOSE","tools","R.utils")

loadpackages = function(packages){
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
}

options(stringsAsFactors = FALSE)

init_params = list()

###### Cluster commands ######
if(!interactive()){
  invisible(suppressMessages(if(!require("R.utils",character.only = TRUE,quietly = TRUE)){
    install.packages("R.utils")
  }))
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  must_args = c("wd","name","anno")
  if(!all(must_args %in% names(args))){
    print_help <- function() {
      title = "ChIP analysis"
      opts = rbind(c("Choose to analyze MACS xls files or raw ChIP FASTQ files"),
                   c("Experimental design file must be Sample name [tab] xls input [tab] xls ChIP (for MACS data)"),
                   c("Experimental design file must be Sample name [tab] FASTQ input [tab] FASTQ ChIP (for FASTQ single-end data)"),
                   c("Experimental design file must be Sample name [tab] FASTQ mate1 input [tab] FASTQ mate2 input [tab] FASTQ mate1 ChIP [tab] FASTQ mate2 ChIP (for FASTQ paired-end data)"),
                   c(""),
                   c("--exp    ","Experimental design file"),
                   c("--wd    ","Working directory path"),
                   c("--name    ","Experiment name"),
                   c("--anno    ","Annotation file path (GFF or GTF file)"),
                   c("--analyze    ","Load an RData file and proceed to analysis"),
                   c("--fastq    ","Performs trimming, alignment and binning of FASTQ data (instead of MACS; default = FALSE)"),
                   c("--genome    ","Genome FASTA file (required for --fastq)"),
                   c("--fraglen    ","Fragment length (for --fastq; Default = 200)"),
                   c("--bin    ","Bin size (for --fastq; Default = 200)"),
                   c("--org    ","Organism number (if exists)"),
                   c("--key    ","Gene key for organism"),
                   c("--t    ","Number of compute threads (Default = detected cores)"),
                   c("--list_orgs    ","List org.db and exit"),
                   c("--arg    ","Additional R arguments (multiple arguments in separate flags)"))
      lines = c("Usage: Rscript ChIPSeq.R [options]","",title,"","Options:",apply(opts, 1, function(r) sprintf("  %-*s  %s", max(nchar(opts[,1])), r[1], r[2]))    )
      cat(paste0(lines, collapse = "\n"), "\n")
    }
    print_help()
    stop(paste0("Missing command line input --> ",paste(must_args[!must_args %in% names(args)],collapse = " | ")), call. = TRUE)
  }
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  
  loadpackages(packages = packages)
  if("list_orgs" %in% names(args)){
    org_list = suppressMessages(BiocManager::available(pattern = "org[.](.*)[.]db",include_installed = TRUE))
    cat(paste0(seq_along(org_list),": ",org_list,"\n"))
    stop("Run script with --org and organism number",call. = TRUE)
  }
  
  init_params[["wd"]] = args[["wd"]]
  init_params[["name"]] = args[["name"]]
  init_params[["annotation"]] = normalizePath(args[["anno"]])
  
  if("org" %in% names(args)){
    init_params[["org_db"]] = suppressMessages(BiocManager::available(pattern = "org[.](.*)[.]db",include_installed = TRUE))[as.numeric(args[["org"]])]
    if(!suppressMessages(require(init_params[["org_db"]],character.only = TRUE,quietly = TRUE))){
      BiocManager::install(init_params[["org_db"]],update = FALSE,ask = FALSE)
      library(init_params[["org_db"]],character.only = TRUE)
    }
    org_db_link = get(init_params[["org_db"]])
    if("key" %in% names(args)){
      init_params[["key"]] = keytypes(org_db_link)[as.numeric(args[["key"]])]
    }else{
      key_list = keytypes(org_db_link)
      cat(paste0(seq_along(key_list),": ",key_list,"\n"))
      stop("Run script with --key and key number",call. = TRUE)
    }
  }else{
    init_params[["org_db"]] = ""
  }
  
  if("analyze" %in% names(args)){
    init_params[["fastq"]] = FALSE
    init_params_rem = init_params
    load(init_params[["analyze"]])
    init_params = init_params_rem
    rm(init_params_rem)
    args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
    genome_annotation = txdbmaker::makeTxDbFromGFF(file = init_params[["annotation"]],format = "auto")
    if("org_db" %in% names(init_params)){
      org_db_link = get(init_params[["org_db"]])
    }
  }else{
    init_params[["exp"]] = normalizePath(args[["exp"]])
    
    init_params[["fastq"]] = "fastq" %in% names(args)
    if(init_params[["fastq"]]){
      if("fraglen" %in% names(args)){
        init_params[["fraglen"]] = as.numeric(args[["fraglen"]])
      }else{
        init_params[["fraglen"]] = 200
      }
      if("bin" %in% names(args)){
        init_params[["bin"]] = as.numeric(args[["bin"]])
      }else{
        init_params[["bin"]] = 200
      }
      if("genome" %in% names(args)){
        init_params[["genome"]] = normalizePath(args[["genome"]])
      }else{
        stop("Please provide a genome FASTA file using --genome", call. = TRUE)
      }
    }
  }
  
  cat("Working directory: ",getwd(),"\n", sep = "")
  cat("Experiment name: ",init_params[["name"]],"\n", sep = "")
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
  
}else{
  loadpackages(packages = packages)
  
init_params[["threads"]] = detectCores()
  
###### Set working directory ######
cat("#####   Select working directory   #####\n")
init_params[["wd"]] = selectDirectory(caption = "Choose working directory")

init_params[["name"]] = as.character(readline(prompt = "Select experiment name: "))
init_params[["exp"]] = selectFile(caption = "Select experimental design file: ",path = getwd())
##### Select files and experimental design ######
init_params[["fastq"]] = grepl(pattern = "^y",x = as.character(readline(prompt = "FASTQ input data (Y/N)? ")),ignore.case = TRUE)
if(init_params[["fastq"]]){
  init_params[["fraglen"]] = as.numeric(readline(prompt = "Choose fragment length (Default = 200): "))
  init_params[["bin"]] = as.numeric(readline(prompt = "Choose bin size (Default = 200): "))
  if(is.na(init_params[["fraglen"]])){init_params[["fraglen"]] = 200}
  if(is.na(init_params[["bin"]])){init_params[["bin"]] = 200}
  init_params[["genome"]] = selectFile(caption = "Select genome file: ",path = getwd())
}
init_params[["annotation"]] = selectFile(caption = "Select GFF/GTF file: ",path = getwd())
init_params[["org_db"]] = select.list(choices = suppressMessages(BiocManager::available(pattern = "org[.](.*)[.]db",include_installed = TRUE)),title = "Select annotation database (0 for none)",multiple = FALSE)
if(!suppressMessages(require(init_params[["org_db"]],character.only = TRUE,quietly = TRUE))){
  BiocManager::install(init_params[["org_db"]],update = FALSE,ask = FALSE)
}
org_db_link = get(init_params[["org_db"]])
if(init_params[["org_db"]] != ""){
  init_params[["key"]] = select.list(choices = keytypes(org_db_link),title = "Select annotation key for genes",multiple = FALSE)
}
}

##### Register threads #####
if(.Platform$OS.type == "unix"){
  register(BPPARAM = MulticoreParam(workers = init_params[["threads"]]))
}else{
  register(BPPARAM = SnowParam(workers = init_params[["threads"]]))
}

##### Input prep #####
if(!exists("experimental_design")){
experimental_design = read.table(file = init_params[["exp"]],header = FALSE,sep = "\t")
if(ncol(experimental_design) == 2){
  colnames(experimental_design) = c("Name","Peaks")
}else{
if(ncol(experimental_design) == 3){
  colnames(experimental_design) = c("Name","Input","ChIP")
}else{
  if(ncol(experimental_design) == 5){
  colnames(experimental_design) = c("Name","Input1","Input2","ChIP1","ChIP2")
  }else{
    stop("Invalid experimental design. Please see help for details.",call. = TRUE)
  }
}
}
for(i in 2:ncol(experimental_design)){
  experimental_design[[i]] = normalizePath(experimental_design[[i]])
}
genome_annotation = txdbmaker::makeTxDbFromGFF(file = init_params[["annotation"]],format = "auto")
if(init_params[["fastq"]]){
  genome_chromosomes = gsub(x = gsub(x = grep(pattern = "^>",x = readLines(init_params[["genome"]]),value = TRUE),pattern = " .*$",replacement = ""),pattern = ">",replacement = "")
if(!all(genome_annotation$user_seqlevels %in% genome_chromosomes)){
  cat("Genome sequence chromosomes: ",genome_chromosomes,"\n")
  cat("Genome annotation chromosomes: ",genome_annotation$user_seqlevels,"\n")
  stop(paste0("Some annotations do not have matching chomosome sequences: ",genome_annotation$user_seqlevels[!genome_annotation$user_seqlevels %in% genome_chromosomes]),call. = TRUE)
}
}
}

print(data.frame(Value = sapply(init_params,function(x) paste(x,collapse = ", "))))
setwd(init_params[["wd"]])

##### FASTQ analysis #####
if(init_params[["fastq"]]){
  init_params[["paired"]] = ncol(experimental_design) > 3
  cat("Data type: ",if(init_params[["paired"]]){"Paired-end"}else{"Single-end"},"\n", sep = "")
start_time=Sys.time()

##### Trim fastq, align and QC ######
cl = makeCluster(init_params[["threads"]])
dir.create("Trimmed",showWarnings = FALSE)
dir.create("Aligned",showWarnings = FALSE)
dir.create("Binned",showWarnings = FALSE)
for(i in unlist(experimental_design[,-1],use.names = FALSE)){
  cat("#####   Trimming FASTQ for ",basename(i),"   #####\n", sep = "")
  fastqFilter(fn = i, fout=paste0("Trimmed/",basename(i),"_trimmed.fastq.gz"),truncQ = 2, truncLen = 0, trimLeft = 0, maxN=0, minQ=0, maxEE=Inf,rm.phix=FALSE, n=1e+06, compress = TRUE, verbose=FALSE)
}

for(i in 2:ncol(experimental_design)){
  experimental_design[[i]] = normalizePath(paste0("Trimmed/",sub(x = basename(experimental_design[[i]]),pattern = "$",replacement = "_trimmed.fastq.gz")))
}
if(init_params[["paired"]]){
  tmp1 = cbind(experimental_design[,c("Input1","Input2")],gsub(pattern = "$",replacement = "_Input",experimental_design[["Name"]]))
  colnames(tmp1) = c("FileName1","FileName2","SampleName")
  tmp2 = cbind(experimental_design[,c("ChIP1","ChIP2")],gsub(pattern = "$",replacement = "_ChIP",experimental_design[["Name"]]))
}else{
  tmp1 = cbind(experimental_design[,"Input",drop = FALSE],gsub(pattern = "$",replacement = "_Input",experimental_design[["Name"]]))
  colnames(tmp1) = c("FileName","SampleName")
  tmp2 = cbind(experimental_design[,"ChIP",drop = FALSE],gsub(pattern = "$",replacement = "_ChIP",experimental_design[["Name"]]))
}
colnames(tmp2) = colnames(tmp1)
runfile = rbind(tmp1,tmp2)
rm(tmp1,tmp2)
write.table(x = runfile,file = "runfile.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)

#### TODO bowtie alignment not working
for(i in 1:nrow(runfile)){
  cat("#####   Aligning ",runfile[i,"SampleName"],"   #####\n", sep = "")
  if(!init_params$paired){
    qAlign(sampleFile = "runfile.txt",genome = init_params[["genome"]],clObj = cl,alignmentsDir = "Aligned")
    file.rename(from = paste0("Aligned/",list.files(pattern=paste0("^",runfile[i,"SampleName"],"(.*).bam$"),path = "Aligned")), to = paste0("Aligned/",runfile[i,"SampleName"],".bam"))
    file.rename(from = paste0("Aligned/",list.files(pattern=paste0("^",runfile[i,"SampleName"],"(.*).bam.bai$"),path = "Aligned")), to = paste0("Aligned/",runfile[i,"SampleName"],".bam.bai"))
    file.rename(from = paste0("Aligned/",list.files(pattern=paste0("^",runfile[i,"SampleName"],"(.*).bam.txt$"),path = "Aligned")), to = paste0("Aligned/",runfile[i,"SampleName"],".bam.txt"))
    cat("#####   Constructing ",init_params[["bin"]],"bp bins for ",runfile[i,"SampleName"],"   #####\n", sep = "")
    constructBins(infile = paste0("Aligned/",runfile[i,"SampleName"],".bam"), fileFormat = "bam", outfileLoc = "Binned/",byChr = FALSE, useChrfile = FALSE,chrfile = NULL, excludeChr = NULL, PET = FALSE, fragLen = init_params[["fraglen"]],binSize = init_params[["bin"]],capping = 0)
  }else{
    qAlign(sampleFile = "runfile.txt",genome = init_params[["genome"]],clObj = cl,alignmentsDir = "Aligned",paired = "fr")
    file.rename(from = paste0("Aligned/",list.files(pattern=paste0("^",runfile[i,"SampleName"],"(.*).bam$"),path = "Aligned")), to=paste0("Aligned/",runfile[i,"SampleName"],".bam"))
    file.rename(from = paste0("Aligned/",list.files(pattern=paste0("^",runfile[i,"SampleName"],"(.*).bam.bai$"),path = "Aligned")), to=paste0("Aligned/",runfile[i,"SampleName"],".bam.bai"))
    file.rename(from = paste0("Aligned/",list.files(pattern=paste0("^",runfile[i,"SampleName"],"(.*).bam.txt$"),path = "Aligned")), to=paste0("Aligned/",runfile[i,"SampleName"],".bam.txt"))
    cat("#####   Constructing ",init_params[["bin"]],"bp bins for ",runfile[i,"SampleName"],"   #####\n", sep = "")
    constructBins(infile = paste0("Aligned/",runfile[i,"SampleName"],".bam"), fileFormat = "bam", outfileLoc = "Binned/",byChr = FALSE, useChrfile = FALSE,chrfile = NULL, excludeChr = NULL, PET = TRUE,binSize = init_params[["bin"]],capping = 0)
  }
  cat("#####   Creating QC report for ",runfile[i,"SampleName"],"   #####\n", sep = "")
  qQCReport(input = "Aligned",pdfFilename = paste0(runfile[i,"SampleName"],"_QCReport.pdf"), clObj = cl)
  unlink("runfile.txt")
}
stopCluster(cl)

##### Fitting to model #####
tmp = as.data.frame(matrix(NA,nrow = nrow(experimental_design),ncol = 2,dimnames = list(NULL,c("Name","Peaks"))))
for(i in 1:nrow(experimental_design)){
  cat("#####   Reading input and ChIP bins for ",experimental_design[i,"Name"],"   #####\n", sep = "")
  if (!paired){
    bin = readBins(type = c("input","chip"),fileName = c(paste0("Binned/",experimental_design[i,"Name"],c("_Input","_ChIP"),"_fragL",init_params[["fraglen"]],"_bin",init_params[["bin"]],".txt")),parallel = TRUE,nCore = init_params[["threads"]]-1)
  }else{
    bin = readBins(type = c("input","chip"),fileName = c(paste0("Binned/",experimental_design[i,"Name"],c("_Input","_ChIP"),"_bin",init_params[["bin"]],".txt")),parallel = TRUE,nCore = init_params[["threads"]]-1)
}

cat("#####   Fitting and creating peaks for ",runfile[i,"SampleName"],"   #####\n", sep = "")
fit = mosaicsFit(bin,analysisType="IO",bgEst="rMOM")
peak = mosaicsPeak(fit, signalModel="2S", FDR=0.05, maxgap=200, minsize=50, thres=10)

##### Exporting results #####
cat("#####   Exporting peaks for ",runfile[i,"SampleName"],"   #####\n", sep = "")
BED_file = normalizePath(paste0("Binned/",experimental_design[i,"Name"],"_bin",init_params[["bin"]],"_Sharp-peak.bed"))
export(peak, type = "bed", filename=BED_file)
write.table(read.table(BED_file, header=FALSE)[-1,], file=gsub(pattern = ".bed$",replacement = "_header-eliminated.bed",x = BED_file), sep="\t", col.names=FALSE, row.names=FALSE)
tmp[i,] = c(experimental_design[i,"Name"],BED_file)
}
experimental_design = tmp
rm(tmp)
save.image(paste0(init_params[["name"]],"_aligned.RData"))
cat("#####  FASTQ normalization, mapping and peak detection run finished in ",format(round(Sys.time()-start_time,2),nsmall=2),"   #####\n", sep = "")
}

##### Analysis #####
start_time=Sys.time()

dir.create("Results/Coverage_plot",showWarnings = FALSE,recursive = TRUE)
dir.create("Results/Gene_groups",showWarnings = FALSE,recursive = TRUE)
dir.create("Results/GO_enrichment",showWarnings = FALSE,recursive = TRUE)
dir.create("Results/Vennpie",showWarnings = FALSE,recursive = TRUE)
dir.create("Results/Annotation_overlap",showWarnings = FALSE,recursive = TRUE)
dir.create("Results/Barchart",showWarnings = FALSE,recursive = TRUE)
dir.create("Results/Pie_chart",showWarnings = FALSE,recursive = TRUE)
dir.create("Results/Heatmap",showWarnings = FALSE,recursive = TRUE)
dir.create("Results/Average_profile",showWarnings = FALSE,recursive = TRUE)
dir.create("Results/Gene_set",showWarnings = FALSE,recursive = TRUE)
if(!init_params[["fastq"]] & !"analyze" %in% names(init_params)){
  init_params[["bin"]] = "_MACS"
  for(i in unlist(experimental_design[["Peaks"]],use.names = FALSE)){
    cat("#####   Formatting MACS2 input for ",i,"   #####\n", sep = "")
    write.table(read.table(file = i,header = TRUE,comment.char = "#")[,c("chr","start","end","name","pileup")],file = sub(x = i,pattern = "[.]xls",replacement = ".bed"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,eol = "\n")
  }
  experimental_design$Peaks = sub(x = unlist(experimental_design[["Peaks"]],use.names = FALSE),pattern = "[.]xls",replacement = ".bed")
}

promoters = getPromoters(TxDb=genome_annotation, upstream=3000, downstream=3000)
tagMatrix = list()
peakAnno = list()
gene = list()

##### Read each peak file #####
for(idx in 1:nrow(experimental_design)){
  bed = experimental_design[idx,"Peaks"]
  bed_name = experimental_design[idx,"Name"]
cat("#####   Reading peak file for ",bed_name,"   #####\n", sep = "")
peak_BED = readPeakFile(peakfile = bed)
tagMatrix[[bed_name]] = getTagMatrix(peak_BED, windows=promoters)
# if(!all(genome_annotation$user_seqlevels %in% genome_chromosomes)){
# chromosome_replace = NULL
# for(j in peak_BED@seqnames@values){
#   chromosome_replace = c(chromosome_replace,chromosome_matrix[match(x = as.character(j),table = chromosome_matrix[,2]),1])
# }
# peak_BED@seqnames@values = factor(x = chromosome_replace,levels = chromosome_replace)
# peak_BED@seqinfo@seqnames = as.character(chromosome_replace)
# }

##### General statistics plots #####
cat("#####   Creating coverage plot for ",bed_name,"   #####\n", sep = "")
png(filename = paste0("Results/Coverage_plot/",bed_name,"_bin",init_params[["bin"]],"_coverage_plot.png"),width = 1440,height = 810,units = "px")
  print(covplot(peak_BED, weightCol = "V5", chrs = genome_annotation$user_seqlevels))
while (!is.null(dev.list())){dev.off()}

cat("#####   Creating average profile for ",bed_name,"   #####\n", sep = "")
png(filename = paste0("Results/Average_profile/",bed_name,"_bin",init_params[["bin"]],"_average_profile.png"),width = 1440,height = 810,units = "px")
  print(plotAvgProf2(peak = peak_BED, conf = 0.95, TxDb = genome_annotation, upstream = 3000, downstream = 3000, xlab = "position", ylab = "read count freq", weightCol = "V5"))
while (!is.null(dev.list())){dev.off()}

cat("#####   Creating heatmap for ",bed_name,"   #####\n", sep = "")
png(filename = paste0("Results/Heatmap/",bed_name,"_bin",init_params[["bin"]],"_heatmap.png"),width = 1440,height = 810,units = "px")
  tagHeatmap(tagMatrix[[bed_name]])
while (!is.null(dev.list())){dev.off()}

cat("#####   Annotating peaks for ",bed_name,"   #####\n", sep = "")
peakAnno[[bed_name]] = annotatePeak(peak = peak_BED, tssRegion = c(-3000, 3000), TxDb = genome_annotation)

cat("#####   Creating piechart for ",bed_name,"   #####\n", sep = "")
png(filename = paste0("Results/Pie_chart/",bed_name,"_bin",init_params[["bin"]],"_piechart.png"),width = 1440,height = 810,units = "px")
  plotAnnoPie(peakAnno[[bed_name]])
while (!is.null(dev.list())){dev.off()}

png(filename = paste0("Results/Barchart/",bed_name,"_bin",init_params[["bin"]],"_barchart.png"),width = 1440,height = 810,units = "px")
  plotAnnoBar(peakAnno[[bed_name]])
while (!is.null(dev.list())){dev.off()}

cat("#####   Creating annotation overlap for ",bed_name,"   #####\n", sep = "")
png(filename = paste0("Results/Annotation_overlap/",bed_name,"_bin",init_params[["bin"]],"_annotation_overlap.png"),width = 1440,height = 810,units = "px")
  upsetplot(peakAnno[[bed_name]], vennpie=FALSE)
while (!is.null(dev.list())){dev.off()}

cat("#####   Creating vennpie for ",bed_name,"   #####\n", sep = "")
png(filename = paste0("Results/Vennpie/",bed_name,"_bin",init_params[["bin"]],"_vennpie.png"),width = 1440,height = 810,units = "px")
  vennpie(peakAnno[[bed_name]])
while (!is.null(dev.list())){dev.off()}

##### Get gene set names #####
cat("#####   Extracting genes for ",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"   #####\n", sep = "")
gene[[bed_name]] = basename(seq2gene(peak_BED, tssRegion=c(-1000,1000), flankDistance=3000, TxDb=genome_annotation))
write.table(x = gene[[bed_name]],file = paste0("Results/Gene_set/",bed_name,"_geneset.txt"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\n")

##### Gene groups annotation #####
# TODO change keyType = "TAIR"
if(init_params[["org_db"]] != ""){
cat("#####   Creating gene groups for ",bed_name,"   #####\n", sep = "")
  #,keyType = "TAIR"
gene_group = groupGO(gene = gene[[bed_name]],OrgDb = org_db_link,keyType = init_params[["key"]],ont = "BP",level = 3,readable = TRUE)
png(filename = paste0("Results/Gene_groups/",bed_name,"_bin",init_params[["bin"]],"_gene_groups.png"),width = 1440,height = 810,units = "px")
  print(barplot(gene_group, drop=TRUE))
while (!is.null(dev.list())){dev.off()}

##### Biological process annotation #####
for(ont in c("BP","MF","CC")){
ont_name = c("biological_process","molecular_function","cellular_compartment")[c("BP","MF","CC") %in% ont]
cat("#####   Creating ",ont_name," gene enrichment for ",bed_name,"   #####\n", sep = "")
gene_enrich = enrichGO(gene = gene[[bed_name]],OrgDb = org_db_link,keyType = init_params[["key"]],ont = ont,qvalueCutoff = 0.05)
if(nrow(gene_enrich) > 0){
  png(filename = paste0("Results/Gene_groups/",bed_name,"_bin",init_params[["bin"]],"_",ont_name,".png"),width = 1440,height = 810,units = "px")
    print(dotplot(gene_enrich))
  while (!is.null(dev.list())){dev.off()}
}
}
}
}

##### Comparative results #####
dir.create("Results/Comparative",showWarnings = FALSE)

cat("#####   Creating average profile overlay   #####\n")
png(filename = paste0("Results/Comparative/Average_profile_overlay_bin",init_params[["bin"]],".png"),width = 1440,height = 810,units = "px")
  print(plotAvgProf(tagMatrix, xlim=c(-3000, 3000)))
while (!is.null(dev.list())){dev.off()}

cat("#####   Creating average plot comparison   #####\n")
png(filename = paste0("Results/Comparative/Average_plot_compare_bin_",init_params[["bin"]],".png"),width = 1440,height = 810,units = "px")
  print(plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row"))
while (!is.null(dev.list())){dev.off()}

cat("#####   Creating heatmap(s)   #####\n")
tag_idx = 1
tag_num = 1
while(tag_idx < nrow(experimental_design)){
  tagsub = as.list(tagMatrix[tag_idx:min(tag_idx + 6,nrow(experimental_design))])
  png(filename = paste0("Results/Comparative/Heatmap_bin",init_params[["bin"]],"_",tag_num,".png"),width = 1440,height = 810,units = "px")
  tagHeatmap(tagsub)
  while (!is.null(dev.list())){dev.off()}
  tag_idx = tag_idx + 7
  tag_num = tag_num + 1
}

cat("#####   Creating annotation bar   #####\n")
png(filename = paste0("Results/Comparative/Annotation_bar_compare_bin",init_params[["bin"]],".png"),width = 1440,height = 810,units = "px")
  plotAnnoBar(peakAnno)
while (!is.null(dev.list())){dev.off()}

##### Comparative annotation #####  
genes_annolist = lapply(peakAnno, function(h) as.data.frame(h)$geneId)
for(ont in c("BP","MF","CC")){
ont_name = c("biological_process","molecular_function","cellular_compartment")[c("BP","MF","CC") %in% ont]
cat("#####   Creating GO enrichment comparison for ",ont_name,"   #####\n", sep = "")
compGO = compareCluster(geneClusters = genes_annolist,fun = "enrichGO",OrgDb = org_db_link,keyType = "TAIR",ont = ont,qvalueCutoff = 0.05)
png(filename = paste0("Results/Comparative/Dotplot_bin",init_params[["bin"]],"_",ont_name,".png"),width = 1440,height = 810,units = "px")
  print(dotplot(compGO))
while (!is.null(dev.list())){dev.off()}
}

save.image(paste0(init_params[["name"]],"_final.RData"))
end_time=Sys.time()
cat("#####  Statistical analysis and results output finished in ",format(round(end_time-start_time,2),nsmall=2),"   #####\n", sep = "")
  
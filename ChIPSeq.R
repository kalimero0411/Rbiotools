##### ChIP analysis #####

packages=c("rChoiceDialogs","BiocParallel","parallel","ShortRead","RSQLite","QuasR","BSgenome",
           "Rsamtools","rtracklayer","GenomicFeatures","Hmisc","Gviz","XML","mosaics","ChIPseeker","clusterProfiler",
           "ReactomePA","dada2","RMariaDB","DOSE","tools","R.utils")
invisible(
  suppressMessages(
    sapply(packages,FUN = function(x) {
      # if(!x %in% rownames(installed.packages())){
      #   cat("Installing package: ",x,"\n", sep = "")
      #   BiocManager::install(x,update = FALSE,ask = FALSE)
      # }
      cat("#####   Loading package: ",x,"   #####\n", sep = "")
      library(x,character.only = TRUE)
    })))

options(stringsAsFactors = FALSE)

###### Cluster commands ######
args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
must_args = c("rdata","wd","name","annotation","t")
if(length(args)){
  if(!all(must_args %in% names(args))){
    help = matrix(data = c("--rdata    ","RData file path",
                           "--wd","Working directory path",
                           "--name","Experiment name (all data will output to a directory by that name in the working directory)",
                           "--annotation","Annotation file path (GFF or GTF file)",
                           "--xls","MACS2 Excel output path (if applicable)",
                           "--paired","Designates that the data is paired-end (if applicable); default = FALSE (single-end)",
                           "--t","Number of compute threads",
                           "--arg","Additional R arguments (multiple arguments in separate flags)")
                  ,ncol = 2,byrow = TRUE)
    prmatrix(help,quote = FALSE,rowlab = rep("",nrow(help)),collab = rep("",2))
    stop(paste0("Missing command line input --> ",paste(must_args[!must_args %in% names(args)],collapse = " | ")), call. = TRUE)
  }
  if(all("MACS2" %in% section,!"xls" %in% names(args))){
    stop(paste0("Missing command line input --> xls ; please specify excel files output from MACS2 using --xls"), call. = TRUE)
  }
  cat("Loading RData file: ",args[["rdata"]],"\n", sep = "")
  load(args[["rdata"]])
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  section = section[!section %in% "Setup"]
  wd = args[["wd"]]
  Experiment_name = args[["name"]]
  genome_annotation_link = args[["annotation"]]
  paired = "paired" %in% names(args)
  dir.create(path = paste0(wd,"/",Experiment_name),showWarnings = FALSE)
  setwd(paste0(wd,"/",Experiment_name))
  cat("Working directory: ",getwd(),"\n", sep = "")
  cat("Experiment name: ",Experiment_name,"\n", sep = "")
  cat("Data type: ",if(paired){"Paired-end"}else{"Single-end"},"\n", sep = "")
  if("MACS2" %in% section){
    xls_files = sub(pattern = ".{1,}/",replacement = paste0(args[["xls"]],"/"),x = xls_files)
  }
  BED_files = NULL
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
cat("#####   Select working directory   #####\n")
wd=rchoose.dir(caption = "Choose working directory:")
setwd(wd)

##### Setup experiment name ######
if(!exists("Experiment_name")){Experiment_name = as.character(readline(prompt = "Select experiment name: "))
if(Experiment_name %in% list.files(wd)){
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
}

##### Select files and experimental design ######
section = rselect.list(choices = c("Setup","FASTQ","MACS2","Analysis"),multiple = TRUE,title = "Select sections to run")
if("Setup" %in% section){
xls_files = NULL
BED_files = NULL
if("FASTQ" %in% section){
  paired = rselect.list(choices = c("Paired-end","Single-end"),multiple = FALSE,title = "Paired or Single end?") == "Paired-end"
    }
if("MACS2" %in% section){
  xls_files=rchoose.files(caption = "Select MACS2 peak.xls files",multi = TRUE)
}
genome_source=rselect.list(choices = c("UCSC","Custom"),multiple = FALSE,title = "Select annotation database")
if(genome_source=="UCSC"){
  genome_annotation_link = rselect.list(choices = as.character(ucscGenomes()[,"db"]),multiple = FALSE,title = "Select UCSC database")
  genome_annotation = makeTxDbFromUCSC(genome = genome_annotation_link)
}
if(genome_source=="Custom"){
  genome_annotation_link = rchoose.files(caption = "Select GFF/GTF file: ")
  genome_annotation = makeTxDbFromGFF(file = genome_annotation_link,format = "auto")
}
org_db = rselect.list(choices = BiocManager::available(pattern = "org[.](.*)[.]db",include_installed = TRUE),multiple = FALSE,title = "Select annotation database")
genome = rchoose.files(caption = "Select genome file: ",multi = FALSE)
genome_chromosomes = sub(x = sub(x = grep(pattern = "^>",x = readLines(genome),value = TRUE),pattern = "[ ](.*)$",replacement = ""),pattern = ">",replacement = "")
Chromosomes = rselect.list(choices = genome_annotation$user_seqlevels,multiple = TRUE,title = "Select chromosomes of interest")
if(!all(genome_annotation$user_seqlevels %in% genome_chromosomes)){
  chromosome_matrix=matrix(data = genome_annotation$user_seqlevels,nrow = length(genome_annotation$user_seqlevels),ncol = 2)
  for(i in seq_along(genome_annotation$user_seqlevels)){
    chromosome_matrix[i,2]=rselect.list(choices = c(as.character(genome_chromosomes),"None of the above"),multiple = FALSE,title = paste0("Select correct chromosome for --> ",genome_annotation$user_seqlevels[i]))
  }
}
}

dir.create(Experiment_name,showWarnings = FALSE)
setwd(paste0(wd,"/",Experiment_name))
if("Setup" %in% section){save.image(paste0("./",Experiment_name,".RData"))}
}

##### FASTQ analysis #####
if("FASTQ" %in% section){
  bin_size = as.integer(readline(prompt = "Select bin size (default = 200): "))
  if(!is.numeric(bin_size)){bin_size = 200}
if(!paired){
  fragment_length = as.integer(readline(prompt = "Select fragment length: "))
  experimental_design = data.frame(Input = rchoose.files(caption = "Select input FASTQ files: ",multi = TRUE),
                                   Experiment = NA_character_,
                                   Name = NA_character_,
                                   stringsAsFactors = FALSE)
  for(i in 1:nrow(experimental_design)){
    experimental_design[i,"Experiment"] = rchoose.files(caption = paste0("Select experiment fastq file for ",basename(experimental_design$Input[i])),multi = FALSE)
    experimental_design[i,"Name"] = readline(prompt = "Choose sample name: ")
  }
}else{
  experimental_design = data.frame(Input_1 = rchoose.files(caption = "Select first mate of each input FASTQ file: ",multi = TRUE),
                                   Input_2 = NA_character_,
                                   Experiment_1 = NA_character_,
                                   Experiment_2 = NA_character_,
                                   Name = NA_character_,
                                   stringsAsFactors = FALSE)
for(i in 1:nrow(experimental_design)){
  experimental_design[i,"Input_2"] = rchoose.files(caption = paste0("Select second pair input fastq file for ",basename(experimental_design[i,"Input_1"])),multi = FALSE)
  experimental_design[i,"Experiment_1"] = rchoose.files(caption = paste0("Select first pair experiment fastq file for ",basename(experimental_design[i,"Input_1"])),multi = FALSE)
  experimental_design[i,"Experiment_2"] = rchoose.files(caption = paste0("Select second pair experiment fastq file for ",basename(experimental_design[i,"Input_1"])),multi = FALSE)
  experimental_design[i,"Name"] = readline(prompt = paste0("Choose sample name for ",paste(basename(as.character(experimental_design[i,c("Input_1","Input_2","Experiment_1","Experiment_2")])),collapse = " | "),": "))
}
}

start_time=Sys.time()
cat("#####   Creating experimental design files   #####\n")
if(!paired){
  experimental_design = data.frame(FileName = c(experimental_design$Input,experimental_design$Experiment),SampleName = c(paste0(experimental_design$Name,"_control"),paste0(experimental_design$Name,"_experiment")))
}else{
  experimental_design = data.frame(FileName1 = c(experimental_design$Input_1,experimental_design$Experiment_1),FileName2 = c(experimental_design$Input_2,experimental_design$Experiment_2),SampleName = c(paste0(experimental_design$Name,"_control"),paste0(experimental_design$Name,"_experiment")))
}


##### Trim fastq, align and QC ######
cl = makeCluster(threads)
experimental_design_trimmed = experimental_design
dir.create("./Trimmed",showWarnings = FALSE)
dir.create("./Aligned",showWarnings = FALSE)
dir.create("./Binned",showWarnings = FALSE)
for(i in unlist(experimental_design[-ncol(experimental_design)],use.names = FALSE)){
  cat("#####   Trimming FASTQ for ",basename(i),"   #####\n", sep = "")
  fastqFilter(fn = i, fout=paste0("./Trimmed/",basename(i),"_trimmed.fastq.gz"),truncQ = 2, truncLen = 0, trimLeft = 0, maxN=0, minQ=0, maxEE=Inf,rm.phix=FALSE, n=1e+06, compress = TRUE, verbose=FALSE)
}

experimental_design_trimmed = data.frame(sub(x = experimental_design[[-ncol(experimental_design)]],pattern = "$",replacement = "_trimmed.fastq.gz"),SampleName = experimental_design[[ncol(experimental_design)]])
for(i in 1:nrow(experimental_design)){
  cat("#####   Aligning ",experimental_design[i,ncol(experimental_design)],"   #####\n", sep = "")
  write.table(x = experimental_design[i,],file = "./runfile.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
  if(!paired){
    qAlign(sampleFile = "./runfile.txt",genome = genome,clObj = cl,alignmentsDir = "./Aligned")
    file.rename(from = paste0("./Aligned/",list.files(pattern=paste0("^",experimental_design[i,ncol(experimental_design)],"(.*).bam$"),path = "./Aligned")), to = paste0("./Aligned/",experimental_design[i,ncol(experimental_design)],".bam"))
    file.rename(from = paste0("./Aligned/",list.files(pattern=paste0("^",experimental_design[i,ncol(experimental_design)],"(.*).bam.bai$"),path = "./Aligned")), to = paste0("./Aligned/",experimental_design[i,ncol(experimental_design)],".bam.bai"))
    file.rename(from = paste0("./Aligned/",list.files(pattern=paste0("^",experimental_design[i,ncol(experimental_design)],"(.*).bam.txt$"),path = "./Aligned")), to = paste0("./Aligned/",experimental_design[i,ncol(experimental_design)],".bam.txt"))
    cat("#####   Constructing ",bin_size,"bp bins for ",experimental_design[i,ncol(experimental_design)],"   #####\n", sep = "")
    constructBins(infile = paste0("./Aligned/",experimental_design[i,ncol(experimental_design)],".bam"), fileFormat = "bam", outfileLoc = "./Binned/",byChr = FALSE, useChrfile = FALSE,chrfile = NULL, excludeChr = NULL, PET = FALSE, fragLen = fragment_length,binSize = bin_size,capping = 0)
  }else{
    qAlign(sampleFile = "./runfile.txt",genome = genome,clObj = cl,alignmentsDir = "./Aligned",paired = "fr")
    file.rename(from = paste0("./Aligned/",list.files(pattern=paste0("^",experimental_design[i,ncol(experimental_design)],"(.*).bam$"),path = "./Aligned")), to=paste0("./Aligned/",experimental_design[i,ncol(experimental_design)],".bam"))
    file.rename(from = paste0("./Aligned/",list.files(pattern=paste0("^",experimental_design[i,ncol(experimental_design)],"(.*).bam.bai$"),path = "./Aligned")), to=paste0("./Aligned/",experimental_design[i,ncol(experimental_design)],".bam.bai"))
    file.rename(from = paste0("./Aligned/",list.files(pattern=paste0("^",experimental_design[i,ncol(experimental_design)],"(.*).bam.txt$"),path = "./Aligned")), to=paste0("./Aligned/",experimental_design[i,ncol(experimental_design)],".bam.txt"))
    cat("#####   Constructing ",bin_size,"bp bins for ",experimental_design[i,ncol(experimental_design)],"   #####\n", sep = "")
    constructBins(infile = paste0("./Aligned/",experimental_design[i,ncol(experimental_design)],".bam"), fileFormat = "bam", outfileLoc = "./Binned/",byChr = FALSE, useChrfile = FALSE,chrfile = NULL, excludeChr = NULL, PET = TRUE,binSize = bin_size,capping = 0)
  }
  cat("#####   Creating QC report for ",experimental_design[i,ncol(experimental_design)],"   #####\n", sep = "")
  qQCReport(input = "./Aligned",pdfFilename = paste0(experimental_design[i,ncol(experimental_design)],"_QCReport.pdf"), clObj = cl)
  file.remove("./runfile.txt")
}
stopCluster(cl)

##### Fitting to model #####
for(i in 1:nrow(experimental_design)){
  cat("#####   Reading input and ChIP bins for ",experimental_design[i,ncol(experimental_design)],"   #####\n", sep = "")
  if (!paired){
    bin = readBins(type = c("input","chip"),fileName = c(paste0("./Binned/",experimental_design[i,ncol(experimental_design)],"_fragL",fragment_length,"_bin",bin_size,".txt")),parallel = TRUE,nCore = threads-1)
  }else{
    bin = readBins(type = c("input","chip"),fileName = c(paste0("./Binned/",experimental_design[i,ncol(experimental_design)],"_bin",bin_size,".txt")),parallel = TRUE,nCore = threads-1)
}

cat("#####   Fitting and creating peaks for ",experimental_design[i,ncol(experimental_design)],"   #####\n", sep = "")
fit = mosaicsFit(bin,analysisType="IO",bgEst="rMOM")
peak = mosaicsPeak(fit, signalModel="2S", FDR=0.05, maxgap=200, minsize=50, thres=10)

##### Exporting results #####
cat("#####   Exporting peaks for ",experimental_design[i,ncol(experimental_design)],"   #####\n", sep = "")
export(peak, type = "bed", filename=paste0("./Binned/",experimental_design[i,ncol(experimental_design)],"_bin",bin_size,"_Sharp-peak.bed"))
write.table(read.table(paste0("./Binned/",experimental_design[i,ncol(experimental_design)],"_bin",bin_size,"_Sharp-peak.bed"), header=FALSE)[-1,], file=paste0("./Results/",experimental_design[i,ncol(experimental_design)],"_bin",bin_size,"_Sharp-peak_header-eliminated.bed"), sep="\t", col.names=FALSE, row.names=FALSE)
BED_files = c(BED_files,paste0("./Binned/",experimental_design[i,ncol(experimental_design)],"_bin",bin_size,"_Sharp-peak.bed"))
}
save.image(paste0("./",Experiment_name,"_aligned.RData"))
end_time=Sys.time()
cat("#####  FASTQ normalization, mapping and peak detection run finished in ",format(round(end_time-start_time,2),nsmall=2),"   #####\n", sep = "")
}

##### Analysis #####
start_time=Sys.time()

if(!"Setup" %in% section){
if(genome_source == "UCSC"){
  genome_annotation = makeTxDbFromUCSC(genome = genome_annotation_link)
}
if(genome_source == "Custom"){
  genome_annotation = makeTxDbFromGFF(file = genome_annotation_link,format = "auto")
}
}

if("MACS2" %in% section){
  experimental_design = data.frame(xls = sub(x = basename(xls_files),pattern = "_peaks.xls",replacement = ""))
  bin_size="_MACS2"
  for(i in xls_files){
    cat("#####   Formatting MACS2 input for ",sub(x = basename(i),pattern = "_peaks.xls",replacement = ""),"   #####\n", sep = "")
    write.table(read.table(file = i,header = TRUE,comment.char = "#")[,c("chr","start","end","name","pileup")],file = sub(x = i,pattern = "[.]xls",replacement = ".bed"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,eol = "\n")
    BED_files = c(BED_files,sub(x = i,pattern = "[.]xls",replacement = ".bed"))
  }
}

if(!org_db %in% rownames(installed.packages())){
  BiocManager::install(org_db,update = FALSE)
}
lapply(org_db,library,character.only = TRUE)
org_db_link = get(org_db)
dir.create("./Results",showWarnings = FALSE)
promoters = getPromoters(TxDb=genome_annotation, upstream=3000, downstream=3000)
tagMatrix_list = as.list(NULL)
peakAnno_list = as.list(NULL)
gene_list = as.list(NULL)

##### Read each peak file #####
if("Analysis" %in% section){
for(bed in BED_files){
cat("#####   Reading peak file for ",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"   #####\n", sep = "")
peak_BED = readPeakFile(peakfile = bed)
if(!all(genome_annotation$user_seqlevels %in% genome_chromosomes)){
chromosome_replace = NULL
for(j in peak_BED@seqnames@values){
  chromosome_replace = c(chromosome_replace,chromosome_matrix[match(x = as.character(j),table = chromosome_matrix[,2]),1])
}
peak_BED@seqnames@values = factor(x = chromosome_replace,levels = chromosome_replace)
peak_BED@seqinfo@seqnames = as.character(chromosome_replace)
}

##### General statistics plots #####
cat("#####   Creating coverage plot for ",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"   #####\n", sep = "")
png(filename = paste0("./Results/",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"_bin",bin_size,"_coverage_plot.png"),width = 1440,height = 810,units = "px")
  print(covplot(peak_BED, weightCol="V5", chrs=Chromosomes))
while (!is.null(dev.list())){dev.off()}
tagMatrix = getTagMatrix(peak_BED, windows=promoters)
tagMatrix_list[[sub(x = basename(bed),pattern = "_peaks.bed",replacement = "")]] = tagMatrix
cat("#####   Creating average profile for ",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"   #####\n", sep = "")
peak_BED_df = data.frame(seqnames=seqnames(peak_BED), starts=start(peak_BED)-1, ends=end(peak_BED), names=peak_BED$V4, score=peak_BED$V5, strand=strand(peak_BED))
dir.create("./peak_data_frames",showWarnings = FALSE)
write.table(peak_BED_df, file=paste0("./peak_data_frames/",basename(file_path_sans_ext(bed)),"_df.bed"), row.names = FALSE, quote = FALSE, col.names = FALSE, sep="\t")
png(filename = paste0("./Results/",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"_bin",bin_size,"_average_profile.png"),width = 1440,height = 810,units = "px")
  print(plotAvgProf2(peak = paste0("./peak_data_frames/",basename(file_path_sans_ext(bed)),"_df.bed"), TxDb = genome_annotation, upstream = 3000, downstream = 3000, xlab = "position", ylab = "read count freq", weightCol = "V5"))
while (!is.null(dev.list())){dev.off()}
cat("#####   Creating heatmap for ",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"   #####\n", sep = "")
png(filename = paste0("./Results/",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"_bin",bin_size,"_heatmap.png"),width = 1440,height = 810,units = "px")
  tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
while (!is.null(dev.list())){dev.off()}
cat("#####   Annotating peaks for ",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"   #####\n", sep = "")
peakAnno = annotatePeak(peak = peak_BED, tssRegion = c(-3000, 3000), TxDb = genome_annotation)
peakAnno_list[[sub(x = basename(bed),pattern = "_peaks.bed",replacement = "")]] = peakAnno
cat("#####   Creating piechart for ",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"   #####\n", sep = "")
png(filename = paste0("./Results/",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"_bin",bin_size,"_piechart.png"),width = 1440,height = 810,units = "px")
  plotAnnoPie(peakAnno)
while (!is.null(dev.list())){dev.off()}
png(filename = paste0("./Results/",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"_bin",bin_size,"_barchart.png"),width = 1440,height = 810,units = "px")
  plotAnnoBar(peakAnno)
while (!is.null(dev.list())){dev.off()}
cat("#####   Creating annotation overlap for ",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"   #####\n", sep = "")
png(filename = paste0("./Results/",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"_bin",bin_size,"_annotation_overlap.png"),width = 1440,height = 810,units = "px")
  upsetplot(peakAnno, vennpie=FALSE,text.scale=3)
while (!is.null(dev.list())){dev.off()}
cat("#####   Creating vennpie for ",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"   #####\n", sep = "")
png(filename = paste0("./Results/",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"_bin",bin_size,"_vennpie.png"),width = 1440,height = 810,units = "px")
  vennpie(peakAnno)
while (!is.null(dev.list())){dev.off()}

##### Get gene set names #####
cat("#####   Extracting genes for ",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"   #####\n", sep = "")
gene=seq2gene(peak_BED, tssRegion=c(-1000,1000), flankDistance=3000, TxDb=genome_annotation)
gene=sub(x = gene,pattern = "^(.*)/",replacement = "")
gene_list[[sub(x = basename(bed),pattern = "_peaks.bed",replacement = "")]] = gene
write.table(x = gene,file = paste0("./Results/",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"_geneset.txt"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\n")

##### Gene groups annotation #####
cat("#####   Creating gene groups for ",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"   #####\n", sep = "")
gene_group=groupGO(gene = gene,OrgDb = org_db_link,keyType = "TAIR",ont = "BP",level = 3,readable = TRUE)
png(filename = paste0("./Results/",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"_bin",bin_size,"_gene_groups.png"),width = 1440,height = 810,units = "px")
  print(barplot(gene_group, drop=TRUE))
while (!is.null(dev.list())){dev.off()}

##### Biological process annotation #####
for(ont in c("BP","MF","CC")){
cat("#####   Creating ",c("biological_process","molecular_function","cellular_compartment")[c("BP","MF","CC") %in% ont]," gene enrichment for ",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"   #####\n", sep = "")
gene_enrich = enrichGO(gene = gene,OrgDb = org_db_link,keyType = "TAIR",ont = ont,qvalueCutoff = 0.05)
if(!is.null(gene_enrich)){
  png(filename = paste0("./Results/",sub(x = basename(bed),pattern = "_peaks.bed",replacement = ""),"_bin",bin_size,"_",c("biological_process","molecular_function","cellular_compartment")[c("BP","MF","CC") %in% ont],".png"),width = 1440,height = 810,units = "px")
    print(dotplot(gene_enrich))
  while (!is.null(dev.list())){dev.off()}
}
}
}

##### Comparative results #####
dir.create("./Results/Comparative",showWarnings = FALSE)
cat("#####   Creating average profile overlay   #####\n")
png(filename = paste0("./Results/Comparative/bin",bin_size,"_average_profile_overlay.png"),width = 1440,height = 810,units = "px")
  print(plotAvgProf(tagMatrix_list, xlim=c(-3000, 3000)))
while (!is.null(dev.list())){dev.off()}
cat("#####   Creating average plot comparison   #####\n")
png(filename = paste0("./Results/Comparative/bin",bin_size,"_average_plot_compare.png"),width = 1440,height = 810,units = "px")
  print(plotAvgProf(tagMatrix_list, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row"))
while (!is.null(dev.list())){dev.off()}
cat("#####   Creating heatmap(s)   #####\n")
for(i in 1:ceiling(nrow(experimental_design)/6)){
  tagsub = as.list(tagMatrix_list[(6*(i-1)+1):min(6*i,nrow(experimental_design))])
  png(filename = paste0("./Results/Comparative/bin",bin_size,"_heatmap_",i,".png"),width = 1440,height = 810,units = "px")
  tagHeatmap(tagsub, xlim=c(-3000, 3000), color=NULL)
  while (!is.null(dev.list())){dev.off()}
}

cat("#####   Creating annotation bar   #####\n")
png(filename = paste0("./Results/Comparative/bin",bin_size,"_annotation_bar_compare.png"),width = 1440,height = 810,units = "px")
  plotAnnoBar(peakAnno_list)
while (!is.null(dev.list())){dev.off()}

##### Comparative annotation #####  
genes_annolist = lapply(peakAnno_list, function(h) as.data.frame(h)$geneId)
for(ont in c("BP","MF","CC")){
cat("#####   Creating GO enrichment comparison for ",c("biological_process","molecular_function","cellular_compartment")[c("BP","MF","CC") %in% ont],"   #####\n", sep = "")
compGO = compareCluster(geneClusters = genes_annolist,fun = "enrichGO",OrgDb = org_db_link,keyType = "TAIR",ont = ont,qvalueCutoff = 0.05)
png(filename = paste0("./Results/Comparative/bin",bin_size,"_",c("biological_process","molecular_function","cellular_compartment")[c("BP","MF","CC") %in% ont],"_dotplot_compare.png"),width = 1440,height = 810,units = "px")
  print(dotplot(compGO))
while (!is.null(dev.list())){dev.off()}
}

}
save.image(paste0("./",Experiment_name,"_final.RData"))
end_time=Sys.time()
cat("#####  Statistical analysis and results output finished in ",format(round(end_time-start_time,2),nsmall=2),"   #####\n", sep = "")
  
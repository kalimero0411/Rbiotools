packages = c("dupRadar","R.utils","parallel")

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

threads = detectCores()
args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
must_args = c("bam","gtf")
if(!all(must_args %in% names(args))){
  help = matrix(data = c("--bam    ","Path to duplicate marked / unmarked BAM file",
                         "--gtf    ","GTF file used in the alignment",
                         "--stranded    ","Strandedness of the FASTQ file (0 = unstranded [default]; 1 = stranded; 2 = reverse)",
                         "--verbose    ","Verbose")
                ,ncol = 2,byrow = TRUE)
  prmatrix(help,quote = FALSE,rowlab = rep("",nrow(help)),collab = rep("",2))
stop(paste0("Missing command line input --> ",paste(must_args[!must_args %in% names(args)],collapse = " | ")), call. = TRUE)
}

if("stranded" %in% names(args)){
  stranded = as.numeric(args[["stranded"]])
}else{
  stranded = 0
}

cat("\n##########\n\n")
if(as.numeric(system(command = paste0("samtools view -f 0x400 -@ ",threads," ",args[["bam"]]," | head | wc -l"),intern = TRUE))){
  cat("Duplicates already marked.\n")
  bamDuprm = args[["bam"]]
}else{
  cat("Marking duplicates...\n")
  bamDuprm = markDuplicates(dupremover="bamutil",
                            bam=normalizePath(args[["bam"]]),
                            path=dirname(system(command = "which bam",intern = TRUE)),
                            rminput=FALSE,
                            threads = threads)
}

cat("Analyzing duplication rates...\n")
dm = analyzeDuprates(bam = bamDuprm,
                     gtf = args[["gtf"]],
                     stranded = stranded,
                     paired = as.numeric(system(command = paste0("samtools view -f 0x1 -@ ",threads," ",args[["bam"]]," | head | wc -l"),intern = TRUE)) > 0,
                     threads = threads,
                     verbose = "verbose" %in% names(args))

cat("Creating plot...\n")
png(filename=paste0(gsub(pattern = "[.][^.]+$",replacement = "",x = bamDuprm),".png"),width=1080,height=1080,units="px")
duprateExpDensPlot(DupMat=dm)
while(!is.null(dev.list())) dev.off()

fit = duprateExpFit(DupMat=dm)
cat("Sample\t",args[["bam"]],"\nIntercept\t",fit$intercept,"\nSlope\t",fit$slope,"\n")
cat("\n##########\n\n")
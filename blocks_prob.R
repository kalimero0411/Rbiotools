packages = c("ggplot2", "dplyr","lpSolve")

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

invisible(suppressMessages(if(!require("R.utils",character.only = TRUE,quietly = TRUE)){
  install.packages("R.utils")
}))
args = R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
must_args = c("blocks", "chr_lengths")
if(!all(must_args %in% names(args))){
  print_help <- function() {
    title = "Brassicaceae Block analysis"
    opts = rbind(c("--blocks    ", "Path to tBLASTn output"),
                 c("--blocks    ", "Path to tBLASTn output"),
                 c("--chr_lengths    ", "Path to chr_sizes file"),
                 c("--match_cutoff    ", "Cutoff for contiguous block size"),
                 c("--chr_cutoff    ", "Cutoff for chromosome size")
    )
    lines = c("Usage: Rscript Blocks_prob.R [options]","",title,"","Options:",apply(opts, 1, function(r) sprintf("  %-*s  %s", max(nchar(opts[,1])), r[1], r[2]))    )
    cat(paste0(lines, collapse = "\n"), "\n")
  }
  print_help()
  stop(paste0("Missing command line input --> ", paste(must_args[!must_args %in% names(args)], collapse = " | ")), call. = TRUE)
}

# query = args[["query"]]
# subject = args[["subject"]]
# system(command = paste0("makeblastdb -in ",normalizePath(subject)," -dbtype nucl"))
# outblast = system(command = paste0("tblastn -query ",normalizePath(query)," -db ",normalizePath(subject)," -evalue 1e-5 -max_target_seqs 100 -outfmt \"6 qacc sacc sstart send length sstrand pident qstart qend evalue\" -num_threads 8"),intern = TRUE)
# blocks = t(as.data.frame(sapply(outblast,function(x) strsplit(x,split = "\t"))))
# rownames(blocks) = NULL
# colnames(blocks) = c("query","subject","sstart","send","length","strand","pident","qstart","qend","evalue")
loadpackages(packages = packages)

blocks = read.table(file = args[["blocks"]], header = FALSE, sep = "\t", col.names = c("Block", "query", "chr", "start", "end", "length", "strand", "pident", "block_start", "block_end", "evalue"))
chr_lengths = read.table(file = args[["chr_lengths"]], header = FALSE, col.names = c("chr", "length"))
if("match_cutoff" %in% names(args)){
  threshold = args[["match_cutoff"]]
}else{
  threshold = 1e6
}
if("chr_cutoff" %in% names(args)){
  chr_threshold = as.numeric(args[["chr_cutoff"]])
}else{
  chr_threshold = 0
}


chr_color = data.frame(Block = unique(blocks$Block), color = "")
chr_color[chr_color$Block %in% c("A", "B", "C"), "color"] = "yellow"
chr_color[chr_color$Block %in% c("D", "E"), "color"] = "red"
chr_color[chr_color$Block %in% c("F", "G", "H"), "color"] = "blue"
chr_color[chr_color$Block %in% c("I", "J"), "color"] = "purple"
chr_color[chr_color$Block %in% c("KL1", "KL2", "MN"), "color"] = "orange"
chr_color[chr_color$Block %in% c("O", "P", "Q", "R"), "color"] = "green"
chr_color[chr_color$Block %in% c("S", "T", "U"), "color"] = "cyan"
chr_color[chr_color$Block %in% c("V", "W", "X"), "color"] = "pink"

chr_lengths = chr_lengths[chr_lengths$length >= chr_threshold, , drop = FALSE]
chr_lengths$chr = factor(chr_lengths$chr)
blocks$chr = factor(blocks$chr, levels = levels(chr_lengths$chr))

block_graph = function(Name,chr_order){
blocks$chr = factor(blocks$chr, levels = chr_order)
chr_lengths$chr = factor(chr_lengths$chr, levels = chr_order)
  
# Create a base plot with chromosome lengths
p = ggplot() +
  # Plot chromosomes as black lines
  geom_segment(data = chr_lengths, aes(y = 0, yend = length, x = chr, xend = chr), color = 'black', linewidth = 1) +
  # Reverse the y-axis to have y=0 at the top
  # Move the x-axis labels to the top
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, max(chr_lengths$length), by = 1e7), labels = function(x) x / 1e6) +
  scale_y_reverse(breaks = seq(0, max(chr_lengths$length), by = 1e7), labels = function(x) x / 1e6) +
  labs(y = "Position (Mbp)", x = "Chromosome") +
  theme(axis.text.x = element_text(face = "bold", size = 15), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15, hjust = 0.5),
        axis.ticks.x.top = element_line(),
        axis.ticks.length.x.top = unit(0.2, "cm"),
        axis.text.x.top = element_text(margin = margin(b = 10), hjust = 1))

for(block_var in unique(blocks$Block)){
  df_positions = blocks[blocks$Block %in% block_var, c("chr", "start", "end")]
  df_positions = df_positions %>%
    group_by(chr) %>%
    arrange(start) %>%
    mutate(group = cumsum(c(1, diff(start) > threshold)))
  
  contiguous_regions = df_positions %>%
    group_by(chr, group) %>%
    summarise(start = min(start), end = max(end), block = block_var,.groups = "drop_last")
  contiguous_regions = contiguous_regions[abs(contiguous_regions$start - contiguous_regions$end) > threshold,]
  
  p = p +
    geom_segment(data = contiguous_regions, aes(y = start, yend = end, x = chr, xend = chr), color = chr_color[chr_color == block_var, "color"], linewidth = 2) +
    geom_text(data = contiguous_regions, aes(y = (start + end) / 2, x = factor(chr), label = block), hjust = -1, color = 'black', size = 7) +
    geom_segment(data = contiguous_regions, aes(y = start, yend = start, x = as.numeric(chr) + 0.1, xend = as.numeric(chr) - 0.1), color = "black", linewidth = 1) +
    geom_segment(data = contiguous_regions, aes(y = end, yend = end, x = as.numeric(chr) + 0.1, xend = as.numeric(chr) - 0.1), color = "black", linewidth = 1)
}

# Print the plot
png(filename = paste0(Name, ".ChrPlot.png"), width = 1080, height = 1080, units = "px")
print(p)
while(!is.null(dev.list())) dev.off()
}

chr_blocks = list()
chr_block_lengths = list()
ACK_blocks = list(
  AK1 = c("A","B","C"),
  AK2 = c("D","E"),
  AK3 = c("F","G","H"),
  AK4 = c("I","J"),
  AK5 = c("MN","KL1","KL2"),
  AK6 = c("O","P","Q","R"),
  AK7 = c("S","T","U"),
  AK8 = c("V","W","X")
)

CEK_blocks = list(
  CK1 = c("A","B"),
  CK2 = c("KL1","KL2","C","D","E"),
  CK3 = c("F","G","H"),
  CK4 = c("MN","J"),
  CK5 = c("U"),
  CK6 = c("I","J","O","P","S","R"),
  CK7 = c("V","W","Q","R","X")
)
geneblock_lengths = c(A = 5053452,B = 3376902,C = 2495088,
                      D = 1778043,E = 3930717,
                      "F" = 6200709,G = 216480,H = 1470926,
                      I = 2483241,J = 4283325,
                      KL1 = 851601,KL2 = 1087645,MN = 4824238,
                      O = 1346490,P = 962199,Q = 1368810,R = 5113284,
                      S = 1483853,"T" = 1060113,U = 6368055,
                      V = 1540245,W = 3114318,X = 1731900)

for(block_var in unique(blocks$Block)){
  df_positions = blocks[blocks$Block %in% block_var,c("chr","start","end")]
  df_positions = df_positions %>%
    group_by(chr) %>%
    arrange(start) %>%
    mutate(group = cumsum(c(1, diff(start) > threshold)))
  
  contiguous_regions = df_positions %>%
    group_by(chr, group) %>%
    summarise(start = min(start), end = max(end), block = block_var)
  contiguous_regions = contiguous_regions[abs(contiguous_regions$start - contiguous_regions$end) > threshold,]
  
  for(i in 1:nrow(contiguous_regions)){
    chr_blocks[[as.character(contiguous_regions[i,"chr"][[1]])]] = c(chr_blocks[[as.character(contiguous_regions[i,"chr"][[1]])]],contiguous_regions[i,"block"])
    chr_block_lengths[[as.character(contiguous_regions[i,"chr"][[1]])]] = c(chr_block_lengths[[as.character(contiguous_regions[i,"chr"][[1]])]],abs(contiguous_regions[i,"end"] - contiguous_regions[i,"start"]))
  }
}
chr_blocks = chr_blocks[order(names(chr_blocks))]
chr_block_lengths = chr_block_lengths[order(names(chr_block_lengths))]
for(i in 1:length(chr_blocks)){
  chr_blocks[[i]] = unlist(chr_blocks[[i]])
}
for(i in 1:length(chr_block_lengths)){
  chr_block_lengths[[i]] = unlist(chr_block_lengths[[i]])
}

Ah_mat = matrix(0,nrow = length(as.character(chr_lengths$chr)),ncol = length(unique(blocks$Block)))
colnames(Ah_mat) = unique(blocks$Block)
rownames(Ah_mat) = as.character(chr_lengths$chr)

for(chr in rownames(Ah_mat)){
  for(block in colnames(Ah_mat)){
    Ah_mat[chr,block] = sum(chr_block_lengths[[chr]][which(chr_blocks[[chr]] %in% block)])
  }
}
Ah_mat = t(apply(Ah_mat,MARGIN = 1,function(x) x/sum(x)))

Chr_prob = function(ancestral_blocks){
  mat = matrix(0,nrow = length(ancestral_blocks),ncol = length(unique(blocks$Block)))
  colnames(mat) = unique(blocks$Block)
  rownames(mat) = names(ancestral_blocks)
  
  for(chr in rownames(mat)){
    for(block in colnames(mat)){
      mat[chr,block] = unname(sum(ancestral_blocks[[chr]] %in% block)*geneblock_lengths[block]/sum(geneblock_lengths[ancestral_blocks[[chr]]]))
    }
  }
  prob_matrix = matrix(0, nrow = nrow(Ah_mat), ncol = nrow(mat))
  for (i in 1:nrow(Ah_mat)) {
    for (j in 1:nrow(mat)) {
      prob_matrix[i, j] = cor(Ah_mat[i, ], mat[j, ])
    }
  }
  prob_matrix = prob_matrix - min(prob_matrix)
  prob_matrix = round(t(apply(prob_matrix, 1, function(x) x / sum(x))),digits = 4)
  colnames(prob_matrix) = rownames(mat)
  rownames(prob_matrix) = rownames(Ah_mat)
  
  prob_vec = as.vector(t(prob_matrix))
  
  # Number of variables
  n_vars1 = nrow(prob_matrix)
  n_vars2 = ncol(prob_matrix)
  
  # Constraints matrix
  f.con = matrix(0, nrow = n_vars1 + n_vars2, ncol = n_vars1 * n_vars2)
  
  # Constraints for AK variables (each can be associated with at most one Chr variable)
  for (i in 1:n_vars1) {
    for (j in 1:n_vars2) {
      f.con[i, (i - 1) * n_vars2 + j] = 1
    }
  }
  
  # Constraints for Chr variables (each must be associated with at least one AK variable)
  for (j in 1:n_vars2) {
    for (i in 1:n_vars1) {
      f.con[n_vars1 + j, (i - 1) * n_vars2 + j] = 1
    }
  }
  
  
  f.rhs = c(rep(1, n_vars1), rep(1, n_vars2))
  f.dir = c(rep("<=", n_vars1), rep(">=", n_vars2))
  result = lp("max", prob_vec, f.con, f.dir, f.rhs, all.bin = TRUE)
  if (result$status == 0) {
    res = data.frame()
    solution = matrix(result$solution, nrow = n_vars1, ncol = n_vars2, byrow = TRUE)
    for (i in 1:n_vars1) {
      for (j in 1:n_vars2) {
        if (solution[i, j] == 1) {
          res = rbind(res,c(rownames(prob_matrix)[i],colnames(prob_matrix)[j],t(prob_matrix)[(i - 1) * n_vars2 + j]))
        }
      }
    }
  } else {
    cat("No optimal solution found.")
  }
  colnames(res) = c("Chromosome","Ancestral","Probability")
  res[["Probability"]] = as.numeric(res[["Probability"]])
  return(res)
}

ACK_prob = Chr_prob(ancestral_blocks = ACK_blocks)
CEK_prob = Chr_prob(ancestral_blocks = CEK_blocks)
write.table(ACK_prob,file = "ACK_probs.table",quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
write.table(CEK_prob,file = "CEK_probs.table",quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

block_graph(Name = args[["blocks"]],chr_order = sort(unique(chr_lengths$chr)))
block_graph(Name = paste0(args[["blocks"]],"_ACK"),chr_order = ACK_prob[order(ACK_prob[["Ancestral"]],-ACK_prob[["Probability"]]),][["Chromosome"]])
block_graph(Name = paste0(args[["blocks"]],"_CEK"),chr_order = CEK_prob[order(CEK_prob[["Ancestral"]],-CEK_prob[["Probability"]]),][["Chromosome"]])
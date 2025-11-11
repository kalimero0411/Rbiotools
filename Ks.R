##### Ks ridge plot from MCScanX output #####
packages = c("ggplot2", "ggridges","tidyverse")

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
must_args = c("gff", "kaks")
if(!all(must_args %in% names(args))){
  print_help <- function() {
    title = "Ks ridge plot from MCScanX output"
    opts = rbind(c("--gff    ", "Path to folder with gff files (.gff / .gff3)"),
               c("--kaks    ", "Path to folder with ka/ks collinearity files (.kaks)"),
               c("--order    ", "Order of species by kaks file names (without extention, comma separated)"),
               c("--color    ", "Species color in order (e.g. 'grey,grey,black,black,black,cyan,cyan')"),
               c("--ks_cutoff    ", "Cutoff for Ks plot (Default = 2)"),
               c("--regex    ","Regex to rename files to samples"))
    lines = c("Usage: Rscript Ks.R [options]","",title,"","Options:",apply(opts, 1, function(r) sprintf("  %-*s  %s", max(nchar(opts[,1])), r[1], r[2]))    )
    cat(paste0(lines, collapse = "\n"), "\n")
  }
  print_help()
  stop(paste0("\nMissing command line input --> ", paste(must_args[!must_args %in% names(args)], collapse = " | ")), call. = TRUE)
}

loadpackages(packages = packages)

if("ks_cutoff" %in% names(args)){
  ks_cutoff = args[["ks_cutoff"]]
}else{
  ks_cutoff = 2
}

if("regex" %in% names(args)){
  name_regex = args[["regex"]]
}else{
  name_regex = ""
}

cat("#####   Loading data   #####\n")
gff_files = list.files(path = args[["gff"]],pattern = "*.gff3?",full.names = TRUE)
gff = lapply(gff_files,function(x){
  read.table(file = x,header = FALSE,sep = "\t",col.names = c("chrom","gene","start","end"))
})
names(gff) = gsub(pattern = ".gff",replacement = "",basename(gff_files))

kaks_files = list.files(path = args[["kaks"]],pattern = "*.kaks$",full.names = TRUE)
kaks_data = lapply(kaks_files,function(x){
  res = read.delim(file = x,header = FALSE,sep = "\t",comment.char = "#")
  colnames(res) = c("block","Gene1","Gene2","pvalue","ka","ks")
  res = res[res[["ka"]] > 0 & res[["ks"]] > 0,,drop = FALSE]
  return(res)
})

names(kaks_data) = gsub(pattern = name_regex,
                        replacement = "",
                        x = gsub(pattern = ".kaks$",replacement = "",basename(kaks_files)))
ks_species = sapply(kaks_data,function(x){
  x[["ks"]]
})

ks_species_data = data.frame(Species = NULL,Ks = NULL)
for(i in 1:length(ks_species)){
  ks_species_data = rbind(ks_species_data,data.frame(Species = names(ks_species)[i],Ks = ks_species[[i]][ks_species[[i]] <= ks_cutoff]))
}

cat("#####   Normalizing values   #####\n")
ks_species_data_normalized = ks_species_data %>%
  group_by(Species) %>%
  mutate(Ks_normalized = Ks) %>%
  ungroup()

ks_species_data_normalized = ks_species_data_normalized %>%
  group_by(Species) %>%
  do({
    density_data = density(.$Ks_normalized, bw = 0.005)
    data.frame(Ks_normalized = density_data$x, Density = density_data$y, Species = .$Species[1])
  }) %>%
  ungroup()

ks_species_data_normalized = ks_species_data_normalized %>%
  group_by(Species) %>%
  mutate(Density_normalized = Density / max(Density)) %>%
  ungroup()

# reorder = c(Outgroup = "Tarenaya_hassleriana",Outgroup = "Carica_papaya",Outgroup = "Capsicum_annuum",
#             Basal = "Aethionema_arabicum",
#             LineageIII = "Euclidium_syriacum",
#             LineageIV = "Arabis_alpina",
#             "Anastatica_hierochuntica",
#             LineageII = "Eutrema_salsugineum",LineageII = "Thlaspi_arvense",LineageII = "Schrenkiella_parvula",LineageII = "Brassica_juncea",LineageII = "Brassica_rapa",LineageII = "Brassica_napus",LineageII = "Brassica_oleracea",
#             LineageI = "Cardamine_hirsuta",LineageI = "Lepidium_sativum",LineageI = "Arabidopsis_thaliana",LineageI = "Arabidopsis_lyrata",LineageI = "Arabidopsis_halleri",LineageI = "Boechera_stricta",LineageI = "Capsella_rubella",LineageI = "Camelina_sativa")

# species_colors = rev(c("grey60","grey60","grey60",
#                    "azure",
#                    "chartreuse2",
#                    "pink",
#                    "white",
#                    "cyan","cyan","cyan","cyan","cyan","cyan","cyan",
#                    "red","red","red","red","red","red","red","red"))
if("order" %in% names(args)){
  reorder = strsplit(args[["order"]],split = ",")[[1]]
}else{
  reorder = sort(unique(ks_species_data_normalized$Species))
}

if("color" %in% names(args)){
  species_colors = rev(strsplit(args[["color"]],split = ",")[[1]])
}else{
  species_colors = rainbow(n = length(unique(ks_species_data_normalized$Species)))
}

# 
# group_colors = c(Outgroup = "grey",
#                  Basal = "purple",
#                  LineageI = "blue",
#                  LineageII = "green",
#                  LineageIII = "red")
# species_colors = sapply(names(reorder),function(x) group_colors[x])
# names(species_colors) = unname(reorder)

cat("#####   Outputting graph   #####\n")
plot = suppressWarnings(ks_species_data_normalized %>%
  mutate(Species = fct_rev(fct_relevel(Species, reorder))) %>%
  ggplot(aes(y = Species, x = Ks_normalized, height = Density_normalized, fill = Species)) +
  geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.6) +
  theme_ridges() +
  theme(
    legend.position = "none",
    panel.spacing = unit(0.01, "lines"),
    strip.text.x = element_text(size = 8),
    axis.title.x = element_text(hjust = 0.5,size = 25),
    axis.title.y = element_text(hjust = 0.5,size = 25),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20,colour = species_colors,hjust = 0)
  ) +
  scale_x_continuous(limits = c(0, ks_cutoff)) +
  scale_fill_manual(values = species_colors) +
  scale_y_discrete(labels = gsub(pattern = "_",replacement = " ",rev(reorder))) +
  xlab("Ks") +
  ylab("Species")
)
png(filename = paste0("Ks_ridgeplot_",format(Sys.time(),'%Y%m%d_%H%M%S'),".png"),width = 1080,height = 1080,units = "px")
print(plot)
while(!is.null(dev.list())) dev.off()

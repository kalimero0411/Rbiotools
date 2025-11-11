packages = c("ggplot2","ggthemes","multcompView","dplyr","R.utils")

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
  must_args = c("files")
  if(!all(must_args %in% names(args))){
    print_help <- function() {
      title = "Create Analysis of Variance plots"
      opts = rbind(c("--files    ","Files of data to analyze"),
                   c("--variable    ","Name for the variable (x axis)"),
                   c("--variable    ","Name for the values (y axis)"))
      lines = c("Usage: Rscript aov_plot.R [options]","",title,"","Options:",apply(opts, 1, function(r) sprintf("  %-*s  %s", max(nchar(opts[,1])), r[1], r[2]))    )
      cat(paste0(lines, collapse = "\n"), "\n")
    }
    print_help()
    stop(paste0("Missing command line input --> ",paste(must_args[!must_args %in% names(args)],collapse = " | ")), call. = TRUE)
  }
  
  loadpackages(packages = packages)
  
  init_params[["files"]] = normalizePath(args[["files"]])
  
  if("variable" %in% names(args)){
    init_params[["variable.name"]] = args[["variable"]]
  }
  
  if("variable" %in% names(args)){
    init_params[["value.name"]] = args[["value"]]
  }
  
}else{
  loadpackages(packages = packages)
  init_params[["files"]] = normalizePath(choose.files(caption = "Choose data files",multi = TRUE))
  init_params[["variable.name"]] = as.character(readline(prompt = "Varaible name: "))
  if(init_params[["variable.name"]] == ""){
    init_params[["variable.name"]] = NULL
  }
  init_params[["value.name"]] = as.character(readline(prompt = "Value name: "))
  if(init_params[["value.name"]] == ""){
    init_params[["value.name"]] = NULL
  }
}

aov_plot = function(input_list,variable.name,value.name){
if(length(input_list) > 2){
  input = stack(lapply(as.list(input_list),function(x) x[!is.na(x)]))
  if(is.null(variable.name)){
    variable.name = "Variable"
  }
  if(is.null(value.name)){
    value.name = "Value"
  }
}else{
  if(all(sapply(input_list,is.numeric) == c(FALSE,TRUE))){
    input = input_list[c(2,1)]
  }else{
    input = input_list
  }
  if(is.null(variable.name)){
    variable.name = names(input_list)[2]
  }
  if(is.null(value.name)){
    value.name = names(input_list)[1]
  }
}
colnames(input) = c("Value","Variable")
anova = aov(Value ~ Variable, data = input)
tukey = TukeyHSD(anova)
dt <- group_by(input, Variable) %>%
  summarise(mean = mean(Value), sd = sd(Value)) %>%
  arrange(desc(mean))
dt$cld <- as.data.frame.list(multcompLetters4(anova, tukey)[[1]])$Letters
val_range = c(min(0,min(dt$mean) - max(dt$sd)),max(dt$mean) + max(dt$sd))

p = ggplot(dt, aes(Variable, mean)) + 
  geom_bar(stat = "identity", aes(fill = mean), show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  labs(x = variable.name, y = paste0("Average ",value.name)) +
  geom_text(aes(label = cld, y = mean + sd + val_range[2]/10)) +
  ylim(val_range*c(1,1.1)) +
  theme_few()

png(filename = paste0(input_file,"_aov.png"),width = 1080,height = 1080,units = "px",res = 150)
print(p)
while(!is.null(dev.list())) dev.off()
}

for(input_file in init_params[["files"]]){
  cat("processing file ",input_file,"...\n")
  input_list = read.table(file = input_file,sep = "\t",header = FALSE)
  aov_plot(input_list = input_list,
           variable.name = init_params[["variable.name"]],
           value.name = init_params[["value.name"]])
}

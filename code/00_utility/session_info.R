# -----------------------------------------------------------------------------------------------------------
# For reproducible research, please install the following R packages 
# and make sure the R and BiocManager versions are correct
# Session Info ----------------------------------------------------------------------------------------------
# setting  value
# version  R version 4.3.1 (2023-06-16)
# os       macOS 15.0.1
# system   x86_64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/New_York
# date     2024-10-22
# rstudio  2024.04.1+748 Chocolate Cosmos (desktop)
# pandoc   3.1.11 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64/ (via rmarkdown)
# -----------------------------------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version = "3.18")

list.of.packages <- c(
  "bacon",
  "coMethDMR",
  "data.table",
  "devtools",
  "DMRcate",
  "doParallel",
  "dorothea",
  "EpiDISH",
  "ExperimentHub", 
  "foreach",
  "GEOquery",                                     
  "ggpubr",  
  "glmnet",
  "grid",
  "gridExtra",
  "GWASTools",  
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  "janitor",
  "limma",
  "lme4",
  "MASS",
  "matrixStats",
  "methylGSA",
  "MethReg",
  "minfi",
  "msigdbr",
  "plyr",   
  "readxl",
  "rGREAT",
  "RPMM",
  "sesame",
  "sesameData",
  "stats",                                        
  "SummarizedExperiment",    
  "survival",
  "survminer",
  "S4Vectors",
  "tidyverse",   
  "wateRmelon",
  "writexl"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  for(new in new.packages){
    if(new %in% available.packages()[,1]){
      install.packages(new)
    } else BiocManager::install(new)
  }
} 


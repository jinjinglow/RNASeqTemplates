# Setup
# ------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2") #Do note that DESeq2 may take a while to run if this is the first time
# install.packages("ggplot2")
# install.packages("matrixStats")
library(DESeq2) 
# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2") #Do note that DESeq2 may take a while to run if this is the first time
# install.packages("ggplot2")
# install.packages("matrixStats")
# install.packages("pheatmap")
library(DESeq2) 
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# If there are any issues with the existing packages, uninstall and install them again
#remove.packages("DESeq2")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2") #Do note that DESeq2 may take a while to run if this is the first time
# ------------------------------------------------------------------------------
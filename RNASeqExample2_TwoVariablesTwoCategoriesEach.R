# ------------------------------------------------------------------------------
# Other examples of design formula:
# - https://www.ebi.ac.uk/gxa/experiments/E-MTAB-8031/Results
#     * RNA-seq on lung tumors
#     * Variables:
#         1. Location: tumor cells or normal tissue nearby
#         2. Sex of the patient
#     * Formula = ~sex + location
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# libraries required
# ------------------------------------------------------------------------------
library(DESeq2) 
library("ggplot2")
library("vsn")
library("pheatmap")
library(EnhancedVolcano)

# ------------------------------------------------------------------------------
# Download data
# ------------------------------------------------------------------------------
library(DESeq2) 
# Download raw counts
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-8031/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
head(counts)
dim(counts) # (61860,27) meaning 61860 genes and 25 samples (2 other columns are geneID and geneName)
# Download metadata
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-8031/resources/ExperimentDesignFile.RnaSeq/experiment-design")
head(metadata)
dim(metadata) # (total 25 rows) 18 different information about the 25 samples

# ------------------------------------------------------------------------------
# Wrangle data for DESeq2
# ------------------------------------------------------------------------------

# DESeq expects the counts to have gene IDs as row names
head(counts)
rownames(counts) = counts$Gene.ID
head(counts)

# Remove unused columns (gene ID and gene name)
genes = counts[, c("Gene.ID", "Gene.Name")]
head(genes) # extract genes columns
counts = counts[, -c(1, 2)] # remove genes columns in counts
head(counts)

# DESeq expects the metadata matrix to have sample IDs in the rownames
head(metadata)
rownames(metadata) = metadata$Run
head(metadata)

# Only keep columns of interest
colnames(metadata) # find the column name of interest
unique(metadata$Sample.Characteristic.sampling.site.)
unique(metadata$Sample.Characteristic.sex.)
metadata = metadata[, c("Sample.Characteristic.sampling.site.", "Sample.Characteristic.sex."), drop = FALSE]

# Look at metadata to see how the variables change with respect to each other
metadata

# Rename column
colnames(metadata) = c("location","sex")
head(metadata)

# Remove spaces in names to avoid DESeq warnings
unique(metadata$location) # print the unique gene names before removal of whitespace
metadata$location <- gsub("\\s+", "", metadata$location)  # removes all kinds of whitespace
unique(metadata$location) # print the unique gene names after removal of whitespace

# ------------------------------------------------------------------------------
# Spot check expression for selected genes
# ------------------------------------------------------------------------------

head(genes) # from this, select a geneName that you wish to draw box plot on
gene_id = genes$Gene.ID[ genes$Gene.Name == 'HRH2' ]
gene_counts = counts[gene_id, ]
gene_counts

gene_data = cbind(metadata, counts=as.numeric(gene_counts))
colnames(gene_data) = c("location","sex","counts")
head(gene_data)

table(gene_data$location) # find out the number of samples
table(gene_data$sex) # find out the number of samples

# Define reference level (or control group)
metadata$sex <- relevel(metadata$sex, ref = "male")
metadata$location <- relevel(metadata$location, ref = "normaltissueadjacenttotumor")
metadata$sex
metadata$location

# ------------------------------------------------------------------------------
# Run DESeq
# ------------------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~ location + sex)


# Ignore genes with low counts
# rowSums(counts(dds) >= 1) >= n → gene must have at least 1 count in n samples
# To Keeps genes that have at least 10 counts in 5 or more samples.
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]

# Run DESeq
dds <- DESeq(dds)

# Extract results for location contrast ("normaltissueadjacenttotumor", "neoplasm")
res <- results(dds, contrast = c("location", "normaltissueadjacenttotumor", "neoplasm"))
# Upregulated genes (log2FC ≥ 1), this means that there is a 2^1=2 times increase in expression
res_upregulated_genes <- subset(res, padj < 0.05 & log2FoldChange >= 1)
# Downregulated genes (log2FC ≤ -1), this means that there is a 2^-1=0.5 times decrease in expression
res_downregulated_genes <- subset(res, padj < 0.05 & log2FoldChange <= -1)
res_upregulated_genes
res_downregulated_genes

cat("Upregulated genes:", nrow(res_upregulated_genes), "\n")
cat("Downregulated genes:", nrow(res_downregulated_genes), "\n")

# Extract results for sex contrast (male vs female)

res <- results(dds, contrast = c("sex", "male", "female"))
# Upregulated genes (log2FC ≥ 1), this means that there is a 2^1=2 times increase in expression
res_upregulated_genes <- subset(res, padj < 0.05 & log2FoldChange >= 1)
# Downregulated genes (log2FC ≤ -1), this means that there is a 2^-1=0.5 times decrease in expression
res_downregulated_genes <- subset(res, padj < 0.05 & log2FoldChange <= -1)
res_upregulated_genes
res_downregulated_genes

cat("Upregulated genes:", nrow(res_upregulated_genes), "\n")
cat("Downregulated genes:", nrow(res_downregulated_genes), "\n")

#-------------------------------------------------------------------------------

# Extract results for sex contrast (male vs female)

res <- results(dds, contrast = c("sex", "male", "female"))
# Upregulated genes (log2FC ≥ 1), this means that there is a 2^1=2 times increase in expression
res_upregulated_genes <- subset(res, padj < 0.05 & log2FoldChange >= 1)
# Downregulated genes (log2FC ≤ -1), this means that there is a 2^-1=0.5 times decrease in expression
res_downregulated_genes <- subset(res, padj < 0.05 & log2FoldChange <= -1)
res_upregulated_genes
res_downregulated_genes

cat("Upregulated genes:", nrow(res_upregulated_genes), "\n")
cat("Downregulated genes:", nrow(res_downregulated_genes), "\n")

# ------------------------------------------------------------------------------
# Plot Results
# ------------------------------------------------------------------------------
# MA plot

#Option 1

# plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples
# Points will be colored blue if the adjusted p value is less than 0.1. 
# Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(res, ylim=c(-2,2))

# Plot counts
# You may specify the specific gene
d <- plotCounts(dds, gene="ENSG00000050555", intgroup="location", 
                returnData=TRUE)

# You may specify the gene index
# d <- plotCounts(dds, gene=which.min(res$padj), intgroup="location", 
#                 returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=location, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


# Option 2

# If you want to change the default significance value for the MA plot, you may use the formula here
# Coerce to a data frame
deseq2ResDF <- as.data.frame(res)
# Examine this data frame
head(deseq2ResDF)
# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA) # choose the level of signicance here
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3)) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Transformation of variance
# ------------------------------------------------------------------------------

# BiocManager::install("vsn")
# this gives log2(n + 1)
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)

#rlog() may take a long time with 50 or more samples, vst() is a much faster transformation
#rld <- rlog(dds, blind=FALSE) 
#meanSdPlot(assay(rld))

library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))

# The plots should look like straight lines,
# which indicates variance is roughly constant across the range of expression 
# If the red line is upward sloping, it is heteroskedastic, not ideal for PCA/clustering

# Heat Map
# install.packages("pheatmap")
library("pheatmap")
ntd <- normTransform(dds) # this gives log2(n + 1)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("location","sex")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Volcano Plots
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')



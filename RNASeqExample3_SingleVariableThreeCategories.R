# ------------------------------------------------------------------------------  
# - https://www.ebi.ac.uk/gxa/experiments/E-GEOD-50760/Results
#     * RNA-seq on human colorectal cancer patients normal colon vs primary tumour vs liver metastases
#     * Variables:
#         1. biopsySite: normal colon vs primary tumour vs liver metastases
#     * Formula = ~biopsySite doesn't make sense since want to compare each genotype separately
#     * Formula = ~group
#         metadata$group = factor(paste(metadata$biopsySite))
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
# From the 'All raw counts for the experiments' in the downloads tab of the above webpage
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-50760/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
head(counts)
dim(counts) # (58735, 56) meaning 58735 genes and 54 samples (2 other columns are geneID and geneName)
# Download metadata
# From the 'Experiment Design' in the downloads tab of the above webpage
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-50760/resources/ExperimentDesignFile.RnaSeq/experiment-design")
tail(metadata)
dim(metadata) # (total 54 rows) 16 different information about the 54 samples

# ------------------------------------------------------------------------------
# Wrangle data for DESeq2
# ------------------------------------------------------------------------------

# DESeq expects the counts to have gene IDs as row names
head(counts)
rownames(counts) = counts$Gene.ID

# Remove unused columns (gene ID and gene name)
genes = counts[, c("Gene.ID", "Gene.Name")]
head(genes)

counts = counts[, -c(1, 2)]
head(counts)

# ------------------------------------------------------------------------------
# Potential additional Data processing before DESeq2
# ------------------------------------------------------------------------------

# 1. Ensure that metadata information only contains records found in counts 
setdiff(metadata$Run,colnames(counts)) # to find which genes in the metadata are missing in counts
# Define sample names from counts
sample_names <- colnames(counts)
# Keep only rows in metadata where Run is in sample_names
metadata <- metadata[metadata$Run %in% sample_names, ]
dim(metadata) # check if there are indeed reduction in rows

# 2. Ensure that the sample names in metadata and counts are arranged properly
# Thus, filter and reorder metadata to match counts
# Else, you will encounter: 'rownames of the colData:are not in the same order as the colnames of the countData:" error
sample_names <- colnames(counts)
metadata <- metadata[match(sample_names, metadata$Run), ]
dim(metadata)

# DESeq expects the metadata matrix to have sample IDs in the rownames
rownames(metadata) <- metadata$Run

# Only keep columns of interest
colnames(metadata) # find the column name of interest
unique(metadata$Factor.Value.biopsy.site.) # there are three sites in total
metadata = metadata[, c("Factor.Value.biopsy.site."), drop = FALSE]

# Rename column
colnames(metadata) = c("biopsySite")
head(metadata)

metadata$biopsySite[metadata$biopsySite == 'primary tumor'] = 'primaryTumor'
metadata$biopsySite[metadata$biopsySite == 'colorectal cancer metastatic in the liver'] = 'liverMetastasis'
unique(metadata$biopsySite)
table(metadata$biopsySite)

# ------------------------------------------------------------------------------
# Spot check expression for selected genes
# ------------------------------------------------------------------------------

head(genes)
gene_id = genes$Gene.ID[ genes$Gene.Name == 'COL11A1' ]
gene_counts = counts[gene_id, ]
head(gene_counts[1:10]) #since there are too many samples, to review the first 10 samples only

gene_data = cbind(metadata, counts=as.numeric(gene_counts))
gene_data

library(ggplot2)
ggplot(gene_data, aes(x = biopsySite, y = counts, fill = biopsySite)) + geom_boxplot() # change to variable name

# Turn genotype into a factor
unique(metadata$biopsySite)
metadata$biopsySite <- factor(metadata$biopsySite, levels = c("normal", "primaryTumor", "liverMetastasis"))
metadata$biopsySite

# Define reference level (or control group)
metadata$biopsySite <- relevel(metadata$biopsySite, ref = "normal")

# ------------------------------------------------------------------------------
# Run DESeq
# ------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~biopsySite)


# Ignore genes with low counts
# rowSums(counts(dds) >= 1) >= n → gene must have at least 1 count in n samples
# To Keeps genes that have at least 10 counts in 5 or more samples.
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]

# Run DESeq
dds <- DESeq(dds)

# Each result from results() gives you one pairwise comparison against the reference group.
# fitting into: log2(count) ~ condition

# primaryTumor vs control
res <- results(dds, contrast = c("biopsySite", "primaryTumor", "normal"))
# Upregulated genes (log2FC ≥ 1), this means that there is a 2^3=8 times increase in expression
res_upregulated_genes <- subset(res, padj < 0.05 & log2FoldChange >= 3)
# Downregulated genes (log2FC ≤ -1), this means that there is a 2^-2=0.125 times decrease in expression
res_downregulated_genes <- subset(res, padj < 0.05 & log2FoldChange <= -3)
res_upregulated_genes
res_downregulated_genes

cat("Upregulated genes:", nrow(res_upregulated_genes), "\n")
cat("Downregulated genes:", nrow(res_downregulated_genes), "\n")

# liverMetastasis vs control
res <- results(dds, contrast = c("biopsySite", "liverMetastasis", "normal"))
# Upregulated genes (log2FC ≥ 1), this means that there is a 2^3=8 times increase in expression
res_upregulated_genes <- subset(res, padj < 0.05 & log2FoldChange >= 3)
# Downregulated genes (log2FC ≤ -1), this means that there is a 2^-2=0.125 times decrease in expression
res_downregulated_genes <- subset(res, padj < 0.05 & log2FoldChange <= -3)
res_upregulated_genes
res_downregulated_genes

cat("Upregulated genes:", nrow(res_upregulated_genes), "\n")
cat("Downregulated genes:", nrow(res_downregulated_genes), "\n")

# ------------------------------------------------------------------------------
# Plot Results
# ------------------------------------------------------------------------------
# MA plot
# plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples
# Points will be colored blue if the adjusted p value is less than 0.1. 
# Points which fall out of the window are plotted as open triangles pointing either up or down.

# Option 1

plotMA(res, ylim=c(-2,2))

# Option 2

# If you want to change the default significance value for the MA plot, you may use the formula here
# Coerce to a data frame
deseq2ResDF <- as.data.frame(res)
# Examine this data frame
head(deseq2ResDF)
# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA) # choose the level of signicance here
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3)) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

# Option 3
# Plot counts for specific gene
d <- plotCounts(dds, gene="ENSG00000006071", intgroup="biopsySite", 
                returnData=TRUE)
# You may specify the gene index
# d <- plotCounts(dds, gene=which.min(res$padj), intgroup="biopsySite", 
#                 returnData=TRUE)

library("ggplot2")
ggplot(d, aes(x=biopsySite, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

# ------------------------------------------------------------------------------
# Transformation of variance
# ------------------------------------------------------------------------------

# BiocManager::install("vsn")
# this gives log2(n + 1)
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))

# rlog() may take a long time with 50 or more samples, vst() is a much faster transformation
# rld <- rlog(dds, blind=FALSE)
# meanSdPlot(assay(rld))

# The plots should look like straight lines,
# which indicates variance is roughly constant across the range of expression 
# If the red line is upward sloping, it is heteroskedastic, not ideal for PCA/clustering

# Heat Map
library("pheatmap")

# Option 1
# Heatmap of top DE genes by ntd
ntd <- normTransform(dds) # this gives log2(n + 1)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20] # define the top n genes
df <- as.data.frame(colData(dds)[,c("biopsySite")])
# Set rownames of df to match column names of the matrix
rownames(df) <- colnames(assay(ntd))
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Option 2
# Heatmap of top DE genes by vsd
topGenes <- head(rownames(res), 50) # define the top n genes
avg_expr <- rowMeans(assay(vsd)[topGenes, ]) # Calculate mean expression per gene
orderedGenes <- topGenes[order(avg_expr, decreasing = TRUE)] #and sort descending
pheatmap(assay(vsd)[orderedGenes,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Volcano Plots
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')


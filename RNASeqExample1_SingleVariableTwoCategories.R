# ------------------------------------------------------------------------------  
# - https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5244/Results
#     * RNA-seq on human breast cancer cells Snai1 knockout vs wildtype
#     * Variables:
#         1. genotype: Snai1 knockout vs wildtype
#     * Formula = ~genotype
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

# Download raw counts
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
head(counts)
dim(counts) # (61860,8) meaning 61860 genes and 6 samples (2 other columns are geneID and geneName)
# Download metadata
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/ExperimentDesignFile.RnaSeq/experiment-design")
head(metadata)
dim(metadata) # (total 6 rows) 18 different information about the 8 samples

# ------------------------------------------------------------------------------
# Wrangle data for DESeq2
# ------------------------------------------------------------------------------

# DESeq expects the counts to have gene IDs as row names
head(counts)
rownames(counts) = counts$Gene.ID
head(counts)

# Remove unused columns (gene ID and gene name)
genes = counts[, c("Gene.ID", "Gene.Name")]
head(genes)
#reason?

counts = counts[, -c(1, 2)]
head(counts)

# DESeq expects the metadata matrix to have sample IDs in the rownames
head(metadata)
rownames(metadata) = metadata$Run
head(metadata)

# Only keep columns of interest
colnames(metadata) # find all the column names and determine the columns of interest
metadata = metadata[, c("Sample.Characteristic.genotype."), drop=FALSE]
# Look at metadata to see how the variables change with respect to each other
metadata

# Rename column
colnames(metadata) = c("genotype") 
metadata

# Remove spaces in names to avoid DESeq warnings
unique(metadata$genotype) # print the unique gene names before removal of whitespace
metadata$genotype <- gsub("\\s+", "", metadata$genotype)  # removes all kinds of whitespace
unique(metadata$genotype) # print the unique gene names after removal of whitespace
# metadata$genotype[metadata$genotype == 'wild type genotype'] = 'wildtype'
# metadata$genotype[metadata$genotype == 'Snai1 knockout'] = 'knockout'
# metadata

# ------------------------------------------------------------------------------
# Spot check expression for knockout gene SNAI1
# ------------------------------------------------------------------------------

head(genes)
# gene_id = genes$Gene.ID[ genes$Gene.Name == 'TSPAN6' ]
# gene_counts = counts[gene_id, ]
# gene_counts
gene_id = genes$Gene.ID[ genes$Gene.Name == 'SNAI1' ]
gene_counts = counts[gene_id, ]
gene_counts

gene_data = cbind(metadata, counts=as.numeric(gene_counts))
gene_data

library(ggplot2)
ggplot(gene_data, aes(x = genotype, y = counts, fill = genotype)) + geom_boxplot()

# Turn genotype into a factor
unique(metadata$genotype) # using the output here, determine your baseline/reference vs the genotype of interest
metadata$genotype = factor(metadata$genotype, levels=c("wildtypegenotype", "Snai1knockout" )) #wildtypegenotype will be baseline
metadata$genotype

# ------------------------------------------------------------------------------
# Run DESeq
# ------------------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~genotype)

# Ignore genes with low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq
dds <- DESeq(dds)

# Compare expression
res = results(dds, contrast=c("genotype", "Snai1knockout", "wildtypegenotype"), alpha=1e-5)
# Upregulated genes (log2FC ≥ 1), this means that there is a 2^1=2 times increase in expression
res_upregulated_genes <- subset(res, padj < 0.05 & log2FoldChange >= 1)
# Downregulated genes (log2FC ≤ -1), this means that there is a 2^-1=0.5 times decrease in expression
res_downregulated_genes <- subset(res, padj < 0.05 & log2FoldChange <= -1)
cat("Upregulated genes:", nrow(res_upregulated_genes), "\n")
cat("Downregulated genes:", nrow(res_downregulated_genes), "\n")
res_upregulated_genes
res_downregulated_genes

# ------------------------------------------------------------------------------
# Plot Results
# ------------------------------------------------------------------------------
# MA plot
# plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples
# Points will be colored blue if the adjusted p value is less than 0.1. 
# Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(res, ylim=c(-2,2))

# Plot counts
# ------------------------------------------------------------------------------
# Option 1: You may specify the specific gene
d <- plotCounts(dds, gene="ENSG00000000003", intgroup="genotype", 
                returnData=TRUE)
# You may specify the gene index
# d <- plotCounts(dds, gene=which.min(res$padj), intgroup="location", 
#                 returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=genotype, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


# Option 2: To display all results, separating significant vs insignificant vy colour
# If you want to change the default significance value for the MA plot, you may use the formula here
# Coerce to a data frame
deseq2ResDF <- as.data.frame(res)
# Examine this data frame
head(deseq2ResDF)
# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA) # choose the level of signicance here
library("ggplot2")
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-3, 3)) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", linewidth=2) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

# ------------------------------------------------------------------------------
# Transformation of variance
# ------------------------------------------------------------------------------

# this gives log2(n + 1)
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))
# The plots should look like straight lines,
# which indicates variance is roughly constant across the range of expression 
# If the red line is upward sloping, it is heteroskedastic, not ideal for PCA/clustering

# Heat Map
library("pheatmap")
ntd <- normTransform(dds) # this gives log2(n + 1)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("genotype")])
# Set rownames of df to match column names of the matrix
rownames(df) <- colnames(assay(ntd))
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

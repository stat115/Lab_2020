#---
# STAT 115 Lab 3 -- 2020
# See longer vignette here:
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
#---

library(stringr)
library(pheatmap)
library(ggplot2)

counts <- readRDS("Lab3_2020_counts.rds")

#-------------
# Process counts into log2tpm
#-------------
log2tpm <- sapply(1:dim(counts)[2], function(idx){
  log2((counts[,idx]/sum(counts[,idx]) * 1000000) + 1)
})
colnames(log2tpm) <- colnames(counts)
sample_data <- data.frame(str_split_fixed(colnames(log2tpm), "_", 3))
colnames(sample_data) <- c("celltype", "replicate", "batch")

# Look for variable genes
rowVars <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
rv <- rowVars(log2tpm)
cutoff <- sort(rv, decreasing = TRUE)[10000]
log2tpm_variable <- log2tpm[rv >= cutoff, ]
pheatmap(cor(log2tpm_variable))

#---------------
# Examine PCA
#---------------
sample_data_pca <- data.frame(
  sample_data,
  prcomp(log2tpm_variable, scale = TRUE, center = TRUE)$rotation
)

ggplot(sample_data_pca, aes(x = PC2, y = PC3, shape = batch, color = celltype)) +
  geom_point(size= 7) 

#---------------
# Now run ComBat
#---------------
library(sva)
combat_dat <- ComBat(log2tpm_variable, sample_data_pca$batch, par.prior = TRUE)
pheatmap(cor(combat_dat))

# Redo PCA
sample_data_pca_combat <- data.frame(
  sample_data,
  prcomp(combat_dat, scale = TRUE, center = TRUE)$rotation
)

ggplot(sample_data_pca_combat, aes(x = PC2, y = PC3, shape = batch, color = celltype)) +
  geom_point(size= 7) 

#---------------
# Now run DESeq2
#---------------
library(DESeq2)
library(apeglm)
ddsMat <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = sample_data,
                                 design = ~ celltype + batch)
ddsMat_wBatch <- DESeqDataSetFromMatrix(countData = counts,
                                         colData = sample_data,
                                         design = ~ celltype )

# Test both adjusting for batch effect and not testing for batch effect.
ddsMat <- DESeq(ddsMat)
ddsMat_wBatch <- DESeq(ddsMat_wBatch)

rdf <- results(ddsMat, contrast = c("celltype", "iPSC", "Fibro"))
rdf_wbatch <- results(ddsMat_wBatch, contrast = c("celltype", "iPSC", "Fibro"))

# Count number of DE genes
sum(rdf$padj < 0.01, na.rm = TRUE)
sum(rdf_wbatch$padj < 0.01, na.rm = TRUE)

# Set up for MA plot
dds_shrink <- lfcShrink(ddsMat, coef="celltype_iPSC_vs_Fibro", type="apeglm")
plotMA(dds_shrink, ylim = c(-5, 5))

# Make data frame for quickly finding the differential genes
library(dplyr)
rdf_df <- data.frame(dds_shrink, gene = rownames(dds_shrink))
rdf_df %>% arrange(padj) %>% head(20)

# Assess direction of the log fold change
log2tpm[c("COL1A2", "LSR"),]


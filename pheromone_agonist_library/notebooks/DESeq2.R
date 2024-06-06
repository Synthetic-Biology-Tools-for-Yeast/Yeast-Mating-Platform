# conda activate yeast_R
# rstudio

# Load necessary libraries
library(tidyverse) ; packageVersion("tidyverse") # 1.3.1
library(phyloseq) ; packageVersion("phyloseq") # 1.22.3
library(vegan) ; packageVersion("vegan") # 2.5.4
library(DESeq2) ; packageVersion("DESeq2") # 1.18.1
library(dendextend) ; packageVersion("dendextend") # 1.10.0
library(viridis) ; packageVersion("viridis") # 0.5.1
library(pheatmap)

setwd("~/Documents/yeast")

# Read the CSV file and set the first column as row names
data <- read.csv('data/interim/finalData/counttable_with_controls.csv', check.names = FALSE, row.names = 1)
sample_info <- read.csv('data/raw/sample_sheet.csv', header=T, row.names = 1, check.names = F)

(condition <- factor(c(rep("FgAgo_gpcr1", 2), rep("FgAgo_gpcr2", 2), rep("FgAgo_gpcr3", 2), rep("FgAgo_gpcr4",2), 
                       rep("BbAgo_gpcr1", 2), rep("BbAgo_gpcr2", 2), rep("BgAgo_gpcr3",2), rep("BbAgo_gpcr4", 2),
                       rep("BcAgo_gpcr1", 2), rep("BcAgo_gpcr2", 2), rep("BcAgo_gpcr3", 2), rep("BcAgo_gpcr4", 2),
                       rep("FoAgo_gpcr1", 2), rep("FoAgo_gpcr2", 2), rep("FoAgo_gpcr3", 2), rep("FoAgo_gpcr4", 2),
                       rep("control_L1", 1), rep("control_L2", 1), rep("control_L3", 1), rep("control_L4", 1))))

(coldata <- data.frame(row.names=colnames(data), condition))

dds <- DESeqDataSetFromMatrix(countData=data, colData=coldata, design=~condition)
nrow(dds)
dds <- dds[ rowSums(counts(dds)) >= 1 ]
nrow(dds)

# Run DESeq
dds <- DESeq(dds)
# Estimate size factors
dds <- estimateSizeFactors(dds)

# Apply variance stabilizing transformation
dds_vst <- varianceStabilizingTransformation(dds)

# Get the transformed count table
vst_trans_count_tab <- assay(dds_vst)

write.csv(vst_trans_count_tab, file = "data/processed/vst_transformed_reads.csv")

# Calculate Euclidian distance matrix
euc_dist <- dist(t(vst_trans_count_tab))

# Perform hierarchical clustering
euc_clust <- hclust(euc_dist, method = "ward.D2")
plot(euc_clust)
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
plot(euc_dend, ylab = "VST Euc. dist.")

# Get differentially abundant features
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Sequence"
head(resdata)
## Write results
write.csv(resdata, file="data/processed/deseq2/diffexpr-results.csv")

# Filter for adjusted p-value < 0.05
sig_res <- res[which(res$padj < 0.05), ]

# Print number of significantly differentially abundant features
print(paste("Number of significantly differentially abundant ASVs: ", nrow(sig_res)))

# Plot heatmap
pheatmap(
  assay(dds)[rownames(sig_res), ], 
  show_rownames = FALSE,
  scale = "row")


########### COMBINED READS ##############
# Read the CSV file and set the first column as row names
data <- read.csv('data/interim/finalData/combined_counts.csv', check.names = FALSE, row.names = 1)

(condition <- factor(c(rep("BbAgo", 2), rep("BcAgo", 2), rep("FgAgo", 2), rep("FoAgo_gpcr4",2), 
                       rep("starting_library", 1))))

(coldata <- data.frame(row.names=colnames(data), condition))

dds <- DESeqDataSetFromMatrix(countData=data, colData=coldata, design=~condition)
nrow(dds)
dds <- dds[ rowSums(counts(dds)) >= 1 ]
nrow(dds)

# Run DESeq
dds <- DESeq(dds)
# Estimate size factors
dds <- estimateSizeFactors(dds)

# Apply variance stabilizing transformation
dds_vst <- varianceStabilizingTransformation(dds)

# Get the transformed count table
vst_trans_count_tab <- assay(dds_vst)

write.csv(vst_trans_count_tab, file = "data/processed/vst_transformed_combined_reads.csv")

# Calculate Euclidian distance matrix
euc_dist <- dist(t(vst_trans_count_tab))

# Perform hierarchical clustering
euc_clust <- hclust(euc_dist, method = "ward.D2")
plot(euc_clust)
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
plot(euc_dend, ylab = "VST Euc. dist.")

# Get differentially abundant features
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Sequence"
head(resdata)
## Write results
write.csv(resdata, file="data/processed/deseq2/diffexpr-results.csv")

# Filter for adjusted p-value < 0.05
sig_res <- res[which(res$padj < 0.05), ]

# Print number of significantly differentially abundant features
print(paste("Number of significantly differentially abundant ASVs: ", nrow(sig_res)))

# Plot heatmap
pheatmap(
  assay(dds)[rownames(sig_res), ], 
  show_rownames = FALSE,
  scale = "row")
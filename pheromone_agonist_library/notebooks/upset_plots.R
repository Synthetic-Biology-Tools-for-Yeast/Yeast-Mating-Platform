# Activate conda environment 'yeast_R'

library(UpSetR)

setwd("~/Documents/yeast")

top5 <- read.csv("data/processed/nt_upset_data.csv", header = TRUE, sep = ",", check.names = F)

# Assuming dataset contains columns '1FgAgo1_A', '1FgAgo2_A', ..., '8FoAgo4_A'
all_sets <- colnames(top5)[2:ncol(top5)]  # Exclude the first column (sequences)

# Open a PDF device with specified width and height
pdf("reports/figures/nt_upset_plot.pdf", width = 15, height = 25)

# Create an UpSet plot with all existing sets
upset(top5, sets = all_sets, sets.bar.color = "#56B4E9", order.by = "freq")

# Close the PDF device
dev.off()



aa_top5 <- read.csv("data/processed/aa_upset_data.csv", header = TRUE, sep = ",", check.names = F)

# Assuming dataset contains columns '1FgAgo1_A', '1FgAgo2_A', ..., '8FoAgo4_A'
all_sets <- colnames(aa_top5)[2:ncol(aa_top5)]  # Exclude the first column (sequences)

# Open a PDF device with specified width and height
pdf("reports/figures/aa_upset_plot.pdf", width = 15, height = 25)

# Create an UpSet plot with all existing sets
upset(aa_top5, sets = all_sets, sets.bar.color = "#56B4E9", order.by = "freq")

# Close the PDF device
dev.off()
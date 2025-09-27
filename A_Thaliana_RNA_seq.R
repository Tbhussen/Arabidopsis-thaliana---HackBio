############################################################
# Differential Expression Analysis in Arabidopsis thaliana
# Using DESeq2
############################################################

# ================================
# 1. Import Libraries
# ================================
suppressPackageStartupMessages({
  library("DESeq2")
  library("pheatmap")
})

# ================================
# 2. Helper Functions
# ================================

# Density plot for p-values
plot_pvalue_density <- function(res, main_title = "P-value Distribution") {
  plot(density(x = na.omit(res$pvalue)),
       main = main_title,
       xlab = "P-value",
       col = "darkgreen",
       lwd = 2)
}

# Volcano plot
plot_volcano <- function(res, padj_cutoff = 0.05, lfc_cutoff = 2, main_title = "Volcano Plot") {
  plot(x = res$log2FoldChange,
       y = -log10(res$padj),
       cex = 0.25,
       pch = 19,
       col = 'grey',
       ylim = c(0, 40),
       ylab = '-log10 Adjusted P-Value',
       xlab = 'Log2 Fold Change',
       main = main_title)
  
  abline(v = c(-lfc_cutoff, lfc_cutoff), h = -log10(padj_cutoff), lwd = 0.5, lty = 2)
  
  upregulated <- subset(res, padj < padj_cutoff & log2FoldChange > lfc_cutoff)
  downregulated <- subset(res, padj < padj_cutoff & log2FoldChange < -lfc_cutoff)
  
  points(upregulated$log2FoldChange,
         -log10(upregulated$padj),
         cex = 0.35, pch = 19, col = 'salmon')
  
  points(downregulated$log2FoldChange,
         -log10(downregulated$padj),
         cex = 0.35, pch = 19, col = 'lightblue')
  
  legend("topright", legend = c("Upregulated", "Downregulated"),
         col = c("salmon", "lightblue"), pch = 19, bty = "n")
  
  return(list(up = upregulated, down = downregulated))
}

# Heatmap of DEGs
plot_heatmap <- function(raw_counts, degs, main_title = "DEG Heatmap") {
  pheatmap(degs,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = FALSE,
           scale = 'row',
           show_colnames = TRUE,
           main = main_title)
}

# ================================
# 3. Set Working Directory
# ================================
setwd("D:/Personal Hub/HackBio")

# ================================
# 4. Load Data
# ================================
A_thaliana_counts <- read.delim("Plant_counts.txt")
A_thaliana_meta   <- read.csv("Plant_metadata.csv")

# Extract raw counts based on metadata samples
raw_counts <- A_thaliana_counts[, as.character(A_thaliana_meta$Sample)]
rownames(raw_counts) <- A_thaliana_counts$Geneid

# ================================
# 5. DESeq2 Analysis
# ================================
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = A_thaliana_meta,
                              design = ~ Type)
dds <- DESeq(dds)
final_res <- results(dds)

# ================================
# 6. Plots
# ================================
plot_pvalue_density(final_res)
res_filtered <- plot_volcano(final_res)

# Heatmap
degs <- rbind(raw_counts[rownames(res_filtered$up), ],
              raw_counts[rownames(res_filtered$down), ])
plot_heatmap(raw_counts, degs)

# ================================
# 7. Export Results
# ================================
write.csv(res_filtered$up, "upregulated_plant.csv")
write.csv(res_filtered$down, "downregulated_plant.csv")
write.csv(raw_counts, "raw_counts_plant.csv")

# Top 100 upregulated
up_top100 <- res_filtered$up[order(res_filtered$up$log2FoldChange, decreasing = TRUE), ][1:100, ]
write.csv(up_top100, "top100_upregulated.csv")

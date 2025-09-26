# Import Libraries
library("DESeq2")
library(pheatmap)

# Set Working Directory
setwd("D:/Personal Hub/HackBio")

# Loading Data
A_thaliana_counts = read.delim("Plant_counts.txt")
A_thaliana_meta = read.csv("Plant_metadata.csv")

head(A_thaliana_counts)
head(A_thaliana_meta)

# Choose Important Columns
raw_counts <- A_thaliana_counts[,as.character(A_thaliana_meta$Sample)]
head(raw_counts)

# Index by gene name
rownames(raw_counts) <- A_thaliana_counts$Geneid
head(raw_counts)

# Create DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = A_thaliana_meta,
                              design = ~Type)
dds
dds$Sample
dds$Type

dds <- DESeq(dds)

final_res <- results(dds)

#we have a truncated data, let's see the distro of p-values
plot(density(x = na.omit(final_res$pvalue)))

#ok let's look at our differentially expressed genes
plot(x = final_res$log2FoldChange, 
     y = -log10(final_res$padj),
     cex = 0.25,
     pch = 19, 
     col = 'grey',
     ylim = c(0,40),
     ylab = 'Adjusted P-Value',
     xlab = 'Log2 FC')

abline(v = c(-2, 2), h = -log10(0.05), lwd = 0.5, lty = 2)
#where are the upregulated
upregulated <- subset(final_res, padj < 0.05 & log2FoldChange > 2)
points(upregulated$log2FoldChange,
       y = -log10(upregulated$padj), 
       cex = 0.35,
       pch = 19,
       col = 'salmon')

#where are the downregulated
downregulated <- subset(final_res, padj < 0.05 & log2FoldChange < -2)
points(downregulated$log2FoldChange,
       y = -log10(downregulated$padj), 
       cex = 0.35,
       pch = 19,
       col = 'lightblue')

mtext('A simple volcano')
#we can merge the two to do a clean and less memory efficient heatmap
degs <- rbind(raw_counts[rownames(upregulated),], 
              raw_counts[rownames(downregulated),])
pheatmap(degs, 
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         scale = 'row',
         show_colnames = T)

#what are the genes that are upregulated
rownames(upregulated)
rownames(downregulated)

#exporting the files
write.csv(upregulated, 'upregulated_plant.csv')
write.csv(downregulated, 'downregulated_plant.csv')
write.csv(raw_counts, 'raw_counts_plant.csv')

up_top100 <- upregulated[order(upregulated$log2FoldChange, decreasing = TRUE), ][1:100, ]
rownames(up_top100)

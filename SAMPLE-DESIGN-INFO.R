count_matrix <- read.csv('E:/BTP-downloads/count_matrix1.csv', row.names = 1)

# Your count matrix
# Your count matrix


# Transpose the count matrix
transposed_matrix <- t(count_matrix)

# Create design matrix
design_matrix <- data.frame(
  Sample_id = colnames(count_matrix),
  Group = c(rep("HF", 8), rep("Control", 11))
)

# Print the design matrix
print(design_matrix)

#
library(DESeq2)
install.packages('pheatmap')
library(pheatmap)
install.packages('RColorBrewer')
library(RColorBrewer)

#set factor level
factors <- factor(design_matrix$Group)
groups <- unique(design_matrix$Group)
groups
groups <- rev(groups)
groups
design_matrix$Group <- factors
design_matrix$Group
#create DSEq object
##count-table=count-matrix and sample-info=design-matrix
dds <- DESeqDataSetFromMatrix(countData = count_matrix,colData = design_matrix,design = ~Group)




#set refernce group

dds$Group <- relevel(dds$Group,ref = "Control")

#genes count atleast 10
keep <- rowSums(counts(dds)>=10)>= min(table(design_matrix$Group))
dds <- dds[keep,]
#perform statistical
#takes a bit time
dds <- DESeq(dds,test = "Wald",sfType = 'poscount')

#get result
deseq_result <- results(dds)

deseq_result
deseq_result <- as.data.frame(deseq_result)
class(deseq_result)



head(deseq_result)


names(deseq_result)



#GENE NAME TO BE HANDELED
# abhi nhi likha hai
deseq_result$Gene-name <- row.names(deseq_result)
head(deseq_result)

#changing some names and subset
deseq_result <- subset(deseq_result,
                       select=c("baseMean","log2FoldChange", "lfcSE","stat","pvalue", "padj" ))

#deseq_result <- subset(deseq_result,select=c("GeneId",and genename,baseMean","log2FoldChange", "lfcSE","stat","pvalue", "padj" ))

names(deseq_result)


#extract de genes with padj <0.05 I changed and log2foldchange <=1 or >=1

deg <- subset(deseq_result,padj>0.05 & abs(log2FoldChange)>=1)


dim(deg)

#92  6

deg <- deg[order(deg$padj),]
head(deg)


#visulasation
#mean normalized count  
#giving error
# Increase the size of the plotting device
# Increase the size of the plotting device
options(repr.plot.width=6, repr.plot.height=6)  # Adjust width and height as needed

# Adjust the margins as needed
par(mar = c(5, 5, 2, 2))  # Format is c(bottom, left, top, right)

# Plot
plotDispEsts(dds, main = "GSE222144 Dispersion Estimates")

# Save the plot to a file
png("dispersion_estimates.png", width = 800, height = 600)  # Adjust width and height as needed
plotDispEsts(dds, main = "GSE222144 Dispersion Estimates")
dev.off()  # Close the png device


#histogram
# Reset the graphics device
dev.off()

# Open a new graphics device with larger dimensions
png("histogram.png", width = 800, height = 600)  # Adjust width and height as needed

# Adjust the margins and reduce font size
par(mar = c(5, 5, 2, 2))  # Adjust margins as needed
hist(deseq_result$padj, breaks = seq(0, 1, length = 11), col = "grey", border = "white", xlab = "", ylab = "", main = "GSE222144 Frequencies of padj-values", cex.axis = 0.8, cex.lab = 0.8)

# Close the graphics device
dev.off()


#volcano

# Reset the graphics device
dev.off()

# Set the palette
old.pal <- palette(c("#00BFFF", "#FF3030"))

# Reset the graphics device
dev.off()

# Adjust the margins and reduce point size
# Open a new graphics device with larger dimensions
png("plot.png", width = 800, height = 600)  # Adjust width and height as needed

# Adjust the margins and reduce point size
par(mar = c(5, 5, 2, 2))  # Adjust margins as needed
plot(deseq_result$log2FoldChange, -log10(deseq_result$padj), 
     main = title, xlab = "log2FC", ylab = "-log10(padj)", pch = 20, cex = 0.3)  # Reduce point size (cex)

# Close the graphics device
dev.off()


#giving colour

#variance
vsd<- vst(dds,blind=FALSE)

#pca plot
plotPCA(vsd,intgroup=c("Group"))

#HEAT MAP of log transformed normalized . we will use top 10 genes



normalized_counts <- counts(dds,normalized=T)

head(normalized_counts)




transfromed_counts <- log2(normalized_counts+1)

head(transfromed_counts)
top_hits <- row.names(deg[1:10,])
top_hits

top_hits<- transfromed_counts[top_hits,]

head(top_hits)
pheatmap(top_hits,cluster_rows = FALSE,cluster_cols = FALSE )



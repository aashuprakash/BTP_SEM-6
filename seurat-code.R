

#"D:\Ayush_Downloads\Gene_GSE235143_merged_counts - _GSE235143_merged_counts.csv"
library(dplyr)
library(Seurat)
library(patchwork)





# Read the gene expression data from the CSV file
pbmc.data <- read.csv("D:/Ayush_Downloads/Gene_GSE235143_merged_counts - _GSE235143_merged_counts.csv", row.names = 1)
#pbmc.data <- read.csv('E:/BTP-downloads/CODE-BATCH/Merged.csv', row.names = 1)

#pbmc.data <- read.csv('E:/BTP-downloads/CODE-BATCH/Merged.csv', row.names = NULL)

# for single dataset
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 5, min.features = 500)
pbmc


# Now the first column name is blank

pbmc.matrix <- as.matrix(pbmc.data)
#pbmc <- CreateSeuratObject(counts = pbmc.matrix, project = "pbmc3k", min.cells = 5, min.features = 300)




# Convert character data to numeric (assuming all values are numeric)
pbmc.data_numeric <- as.matrix(as.data.frame(lapply(pbmc.data, as.numeric)))

# Check for NA values
if (any(is.na(pbmc.data_numeric))) {
  # Handle missing values (if any)
  # For example, you can impute them with mean or median
  pbmc.data_numeric[is.na(pbmc.data_numeric)] <- mean(pbmc.data_numeric, na.rm = TRUE)
}

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data_numeric, project = "pbmc3k", min.cells = 5, min.features = 150)#5,300


# Initialize the Seurat object with the raw (non-normalized) data
#pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 5, min.features = 300)

# Print the summary of the Seurat object
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2



# values changed  wla niche 

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

library(ggplot2)

# Create the VariableFeaturePlot without adjusting size
plot1 <- VariableFeaturePlot(pbmc)

# Add labels
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

# Save plot to an external file (e.g., PNG) with adjusted dimensions
ggsave("variable_features_plot.png", plot1 + plot2, width = 10, height = 8, dpi = 300)



all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#changed values ----

head(pbmc)

pbmc <- subset(pbmc, subset = nFeature_RNA >5000 & nFeature_RNA <20000 & percent.mt >=0)
dim(pbmc)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
dim(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# Adjust plot size
plot1 <- VariableFeaturePlot(pbmc, pt.size = 1)  # Reduce point size to make the plot smaller

# Add labels
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

# Combine plots
combined_plot <- plot1 + plot2

# Display the combined plot
combined_plot


ggsave("combined_plot1.png", combined_plot, width = 10, height = 8)



all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#########GIVING ERROR
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


# Assuming you want to compute 50 principal components
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)


#######

# Check the dimensions of the data matrix
# Check the dimensions of the count data
dim(counts(pbmc))

# Determine the maximum number of principal components
max_pcs <- min(dim(counts(pbmc))) - 1  # Subtract 1 to ensure it's strictly less than min

# Set the number of principal components to a value within the bounds
npcs <- min(max_pcs, 20)  # For example, set to 20 or adjust as needed

# Run PCA with the corrected number of principal components
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = npcs)




#
######error
# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca") + NoLegend()
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)
#############


pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

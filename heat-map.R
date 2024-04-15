

count_data <- read.csv('E:/BTP-downloads/count_matrix1.csv', header=TRUE, row.names=1)


colnames(count_data)

head(count_data)

# Define column names for control and test samples
control_column_names <- paste("HF", 1:8, sep = "_")#HF
test_column_names <- paste("Control", 1:11, sep = "_")

# Rename columns in the count data matrix
colnames(count_data)[1:8] <- control_column_names#HF
colnames(count_data)[9:19] <- test_column_names#Control
count_data <- count_data[which(rowSums(count_data) > 30), ]

count_data
#C-->IS FOR CONTROL AND S--> FOR HF
condtion <- factor(c("HF","HF","HF","HF","HF","HF","HF","HF","Control","Control","Control","Control","Control","Control",'Control',"Control","Control","Control","Control"))


coldata <- data.frame(row.names = colnames(count_data),condtion)
coldata
install.packages("DESeq2")
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~condtion)

dds <- DESeq(dds)

#plots 


vsdata <- vst(dds,blind = FALSE)
#df <- vst(dds,blind = FALSE)
plotPCA(vsdata,intgroup="condtion")
#plotPCA(df,intgroup="condtion")
# Set smaller margins
par(mar = c(5, 4, 4, 2) + 0.1)

# Plot dispersion estimates
plotDispEsts(dds)
res <- results(dds,contrast = c("condtion","Control","HF"))
#head(res)
res
#
df <- res


##########
install.packages("org.Hs.eg.db")

library(org.Hs.eg.db)

res.df < as.data.frame(res)

res.df$symbol <- mapIds(org.Hs.eg.db, keys rownames(res.df), keytype = "ENSEMBL", column = "SYMBOL")
res.df


#####

if(!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
install.packages('texshaping')
library(EnhancedVolcano)
# added parameters 
res  <- na.omit(res)
res <- res[res$baseMean>50,]

EnhancedVolcano(res,x="log2FoldChange",y="padj",lab="ayush")
#
sigs <- na.omit(res)
# adjusting the values
#log2 fold change (MLE): condtion C vs S 
#Wald test p-value: condtion C vs S 
#DataFrame with 0 rows and 6 columns

####value >
#sigs <- sigs[sigs$padj <  0.9999629,]

sigs <- sigs[sigs$padj >  0.05,]


#
sigs
#
df <- as.data.frame(sigs)
df


df.top <- df[ (df$baseMean > 0) & (abs(df$log2FoldChange) > 0),]
df.top

#
rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object
mat<-assay(rlog_out)[rownames(df.top), rownames(coldata)] #sig genes x samples
colnames(mat) <- rownames(coldata)
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mat)


#just kept 25 upper and below 25 values
rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object
mat<-assay(rlog_out)[rownames(df.top), rownames(coldata)] #sig genes x samples
colnames(mat) <- rownames(coldata)
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mat)


##########

nrow(df.top)


#########NEED TO BE CHECKED WHY 2 ?
num_keep <- 70
#1 to num_keep len-num_keep to len
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)) )

#Error in h(simpleError(msg, call)) : 
#error in evaluating the argument 'x' in selecting a method for function 'as.matrix': error in evaluating the argument 'i' in selecting a method for function '[': object 'rows_keep' not found
l2_val <- as.matrix(df.top[rows_keep,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val)<-"logFC"

mean <- as.matrix(df.top[rows_keep,]$baseMean) #getting mean value for each gene we are keeping
colnames(mean)<-"AveExpr"


#DOWNLOADING
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

#
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#
#maps values between b/w/r for min and max l2 values
col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red")) 

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))


#
ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F, 
              column_labels = colnames(mat.scaled), name="Z-score",
              cluster_columns = T)
h2 <- Heatmap(l2_val, row_labels = df.top$symbol[rows_keep], 
              cluster_rows = F, name="logFC", top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(l2_val[i, j],2), x, y)
              })
h3 <- Heatmap(mean, row_labels = df.top$symbol[rows_keep], 
              cluster_rows = F, name = "AveExpr", col=col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(mean[i, j],2), x, y)
              })

h<-h1+h2+h3
h


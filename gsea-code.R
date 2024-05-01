df <- read.csv('E:/BTP-downloads/CODE-BATCH/processed-merged.csv', header=TRUE)

#"E:\BTP-downloads\CODE-BATCH\processed-merged.csv"
# Read in data ===================================================

library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations
install.packages("ggupset")


# Get the genes that are present in your dataframe

genes_in_data <- df$Symbol

#"E:\BTP-downloads\CODE-BATCH\c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt"

#"E:\BTP-downloads\CODE-BATCH\c2.cp.kegg_medicus.v2023.2.Hs.symbols (1).gmt"

#"E:\BTP-downloads\CODE-BATCH\GSEA\c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt"

# Read in the .gmt file
#file <- "E:/BTP-downloads/CODE-BATCH/GSEA/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt"
file <- "E:/BTP-downloads/CODE-BATCH/GSEA/c2.cp.reactome.v2023.2.Hs.symbols.gmt"
pwl2 <- read.gmt(file) 
# Subset to the genes that are present in our dataset
pwl2 <- pwl2[pwl2$gene %in% genes_in_data,] 
# Save the filtered background gene set
filename <- 'reactome.RDS'
saveRDS(pwl2, filename)


# SQUIDTIP! If you want to parse several .gmt files at once, you can use a loop:
#loop things don't run only file <- wla single single run
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
for (file in gmt_files){
  file <- gmt_files[1]
  pwl2 <- read.gmt(file) 
  pwl2 <- pwl2[pwl2$gene %in% genes_in_data,]
  filename <- paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', file)), '.RDS', sep = '')
  saveRDS(pwl2, filename)
}

# Annotate according to differential expression
df <- df %>% mutate(diffexpressed = case_when(
  log2FoldChange > 0 & padj < 0.05 ~ 'UP',
  log2FoldChange < 0 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
))





df <- df[df$diffexpressed != 'NO',]
deg_results_list <- split(df,df$diffexpressed)
# Substitute names so they are annotated nicely in the heatmap later
df$diffexpressed <- gsub('DOWN', 'Healthy', gsub('UP', 'Severe', df$diffexpressed))
unique(df$diffexpressed)
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
deg_results_list <- split(df, df$diffexpressed)


#E:\BTP-downloads\CODE-BATCH\GSEA
## Run ClusterProfiler -----------------------------------------------
bg_path <- "E:/BTP-downloads/CODE-BATCH/GSEA"
"E:/BTP-downloads/CODE-BATCH/GSEA/kegg_legav2023.2.Hs.symbols.gmt.RDS"
# Settings
name_of_comparison <- 'severevshealthy' # for our filename
#background_genes <- 'KEGG' # for our filename
background_genes <- 'reactome'
#bg_genes <- readRDS("E:/BTP-downloads/CODE-BATCH/GSEA/kegg_legav2023.2.Hs.symbols.gmt.RDS") # read in the background genes
bg_genes <- readRDS("E:/BTP-downloads/CODE-BATCH/GSEA/reactome.v2023.2.Hs.symbols.gmt.RDS")
#"E:\BTP-downloads\CODE-BATCH\GSEA\reactome.v2023.2.Hs.symbols.gmt.RDS"
padj_cutoff <- 0.05 # p-adjusted threshold, used to filter out pathways
genecount_cutoff <- 5 # minimum number of genes in the pathway, used to filter out pathways
#filename <- paste0(out_path, 'clusterProfiler/', name_of_comparison, '_', background_genes) # filename of our PEA results
#E:\BTP-downloads\CODE-BATCH\Results-output
out_path <- "E:/BTP-downloads/CODE-BATCH/Results-output"
filename <- paste0(out_path,'clusterProfiler/', name_of_comparison, '_', background_genes)

#not used below

# SQUIDTIP! An option to read in your background genes by only defining your 'background_genes' variable
if(background_genes == 'KEGG'){
  bg_genes <- readRDS(paste0(bg_path, 'kegg.RDS'))
} else if(background_genes == 'reactome'){
  bg_genes <- readRDS(paste0(bg_path, 'reactome.RDS'))
} else if(background_genes == 'go.bp'){
  bg_genes <- readRDS(paste0(bg_path, 'go.bp.RDS'))
} else {
  stop('Invalid background genes. Select one of the following: KEGG, Reactome, GO, or add new pwl to function')
}
#

# Run clusterProfiler on each sub-dataframe
res <- lapply(names(deg_results_list),
              function(x) enricher(gene = deg_results_list[[x]]$Symbol,
                                   TERM2GENE = bg_genes))
names(res) <- names(deg_results_list)

#Convert the list of enrichResults for each sample_pattern to a dataframe with the pathways
res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df) <- names(res)
res_df <- do.call(rbind, res_df)
head(res_df)




#Convert the enrichResults to a dataframe with the pathways
res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df) <- names(res)
res_df <- do.call(rbind, res_df)
head(res_df)


res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust),
                            diffexpressed = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', rownames(res_df)))

# Subset to those pathways that have p adj < cutoff and gene count > cutoff (you can also do this in the enricher function)
target_pws <- unique(res_df$ID[res_df$p.adjust < padj_cutoff & res_df$Count > genecount_cutoff]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
res_df <- res_df[res_df$ID %in% target_pws, ]


print('Saving clusterprofiler results')
#write.csv(res_df, paste0(filename, '_resclusterp.csv'), row.names = FALSE)
# Save the dataframe to a CSV file
write.csv(res_df, "resclusterp_dataframe-reactome.csv", row.names = FALSE)

# VISUALISATION  

#res_df <- read.csv(paste0(out_path, 'clusterProfiler/', 'severevshealthy_reactome_resclusterp.csv'))
#bg_genes <- readRDS(paste0(bg_path, 'reactome.RDS'))

#bg_genes <- readRDS(paste0(bg_path, 'reactome.RDS'))

# Convert clusterProfiler object to a new "enrichResult" object
# Select only upregulated genes in Severe

#"E:\BTP-downloads\R-new-repo\btp1\BTP\resclusterp_dataframe.csv"

originalresdf <- read.csv('E:/BTP-downloads/R-new-repo/btp1/BTP/resclusterp_dataframe.csv', header=TRUE)
````

original_res_df <- res_df

res_df <- res_df %>% filter(diffexpressed == 'UP') %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_df) <- res_df$ID

# For visualisation purposes, let's shorten the pathway names
res_df$Description <- gsub('(H|h)iv', 'HIV', 
                           gsub('pd 1', 'PD-1',
                                gsub('ecm', 'ECM', 
                                     gsub('(I|i)nterleukin', 'IL', 
                                          gsub('(R|r)na', 'RNA', 
                                               gsub('(D|d)na', 'DNA',
                                                    gsub(' i ', ' I ', 
                                                         gsub('(A|a)tp ', 'ATP ', 
                                                              gsub('(N|n)adh ', 'NADH ', 
                                                                   gsub('(N|n)ad ', 'NAD ',
                                                                        gsub('t cell', 'T cell',
                                                                             gsub('b cell', 'B cell',
                                                                                  gsub('built from .*', ' (...)',
                                                                                       gsub('mhc', 'MHC',
                                                                                            gsub('mhc class i', 'MHC I', 
                                                                                                 gsub('mhc class ii', 'MHC II', 
                                                                                                      stringr::str_to_sentence(
                                                                                                        gsub('_', ' ',  
                                                                                                             gsub('GOBP_|KEGG_|REACTOME_', '', res_df$Description)))))))))))))))))))  
``````





originalresdf <- originalresdf %>% filter(diffexpressed == 'UP')  
 # dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(originalresdf) <- originalresdf$ID

# For visualisation purposes, let's shorten the pathway names
original_res_df$Description <- gsub('(H|h)iv', 'HIV', 
                           gsub('pd 1', 'PD-1',
                                gsub('ecm', 'ECM', 
                                     gsub('(I|i)nterleukin', 'IL', 
                                          gsub('(R|r)na', 'RNA', 
                                               gsub('(D|d)na', 'DNA',
                                                    gsub(' i ', ' I ', 
                                                         gsub('(A|a)tp ', 'ATP ', 
                                                              gsub('(N|n)adh ', 'NADH ', 
                                                                   gsub('(N|n)ad ', 'NAD ',
                                                                        gsub('t cell', 'T cell',
                                                                             gsub('b cell', 'B cell',
                                                                                  gsub('built from .*', ' (...)',
                                                                                       gsub('mhc', 'MHC',
                                                                                            gsub('mhc class i', 'MHC I', 
                                                                                                 gsub('mhc class ii', 'MHC II', 
                                                                                                      stringr::str_to_sentence(
                                                                                                        
                                                                                                          gsub('_', ' ',  
                                                                                                               gsub('GOBP_|KEGG_|REACTOME_', '', original_res_df$Description)))))))))))))))))))  

````


enrichres <- new("enrichResult",
                 readable = FALSE,
                 result = original_res_df,################change
                 pvalueCutoff = 0.7,#0.80.05
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.5,#0.4
                 organism = "human",
                 ontology = "UNKNOWN",
                 gene = df$Symbol,
                 keytype = "UNKNOWN",
                 universe = unique(bg_genes$gene),
                 gene2Symbol = character(0),
                 geneSets = bg_genes)
class(enrichres)

# Barplot
barplot(enrichres, showCategory = 20) 
mutate(enrichres, qscore = -log(p.adjust, base = 10)) %>% 
  barplot(x = "qscore")


# Dotplot
dotplot(enrichres, showCategory = 15) + ggtitle("Severe vs Healthy")
# Cnetplot
cnetplot(enrichres)
# Heatplot
heatplot(enrichres, showCategory = 5)
# Treeplot
enrichres2 <- pairwise_termsim(enrichres) # calculate pairwise similarities of the enriched terms using Jaccard’s similarity index
treeplot(enrichres2)

# Enrichment map 
#not giving
emapplot(enrichres2)
# Upsetplot
upsetplot(enrichres)


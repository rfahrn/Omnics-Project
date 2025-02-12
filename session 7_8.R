##################################################################
##################################################################
##################################################################
###                                                            ###
###                                                            ###
###            R Code for Omics for Non-Biologists             ###
###                          Session 7-8                       ###
###                      Date:15 Oct 2024                      ###
###                                                            ###
###                                                            ###
##################################################################
##################################################################
##################################################################


##################################################################
##                     Session Preparation:                     ##
##                    load and attach add-on packages           ##
##################################################################

library(readr)          # fast reading of rectangular data (.csv)
library(DESeq2)         # identification of differentially expressed genes 
library(ggplot2)        # data visualization, complex plots
library(pheatmap)       # data visualization, customizable heatmap
library(EnhancedVolcano)# data visualization, Volcano plots
library(ashr)           # for shrinking estimates, avoids overestimation
library(bannerCommenter)# for commenting 

###===============================================================
###                    Session Preparation:                    ===
###                      Creating Banners                      ===
###===============================================================

banner("R Code for Omics for Non-Biologists \n Session 7-8",
       "Date:15 Oct 2024", 
       emph = F,
       bandChar="#", 
       numLines =3,
       leftSideHashes=3)

banner("Session Preparation:",
       "Creating Banners",
       emph = F,
       bandChar="=", 
       numLines =1,
       leftSideHashes=3)


banner("Session Preparation:",
       "Creating subfolders",
       emph = F,
       bandChar="=", 
       numLines =1,
       leftSideHashes=3)

###===============================================================
###                    Session Preparation:                    ===
###                    Creating subfolders                     ===
###===============================================================
subfolder_names <- c("1_rawdata",
                     "2_metadata",
                     "3_references",
                     "4_results",
                     "5_logs",
                     "6_script",
                     "7_RData")

# Get the current working directory
current_directory <- getwd()

# Create subfolders in the current working directory
for (j in 1:length(subfolder_names)) {
  folder <- dir.create(file.path(current_directory, subfolder_names[j]))
}

# Set the file path
file_path <- "1_rawdata/Mov10_full_counts.txt"

# Read the file with tab delimiter
data <- read_delim(file_path, delim = "\t")

# View the first few rows
head(data)

# ------------------------------------------------------------------------------
# Subset Dataset 

# Subset the data to only the specified columns
subset_data <- data[, c("GeneSymbol", "Mov10_oe_2", "Mov10_oe_3", "Irrel_kd_1", "Irrel_kd_2", "Irrel_kd_3")]


# Prepare the counts data and sample information (metadata)
count_data <- as.matrix(subset_data[, -1]) # exclude GeneSymbol column
rownames(count_data) <- subset_data$GeneSymbol

# Create sample information data frame
sample_info <- data.frame(
 row.names = colnames(count_data),
 condition = c("Mov10_oe", "Mov10_oe", "Irrel_kd", "Irrel_kd", "Irrel_kd")
)

# ------------------------------------------------------------------------------
# Detect Genes using DESEQ2 and Shrink Estimates 

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info, design = ~ condition)

# Pre-filtering: remove rows with very low counts
dds <- dds[rowSums(counts(dds)) > 1, ]

# Run the DESeq pipeline
dds <- DESeq(dds)

# Get the results for Mov10_oe vs Irrel_kd
res <- results(dds, contrast = c("condition", "Mov10_oe", "Irrel_kd"))


# Shrink log fold changes to improve visualization
res_shrink <- lfcShrink(dds, coef = "condition_Mov10_oe_vs_Irrel_kd", type = "ashr")

# ------------------------------------------------------------------------------
# Volcano Plot, PCA plot 

# Create a volcano plot using EnhancedVolcano
EnhancedVolcano(res_shrink,
                lab = rownames(res_shrink),
                x = "log2FoldChange",
                y = "pvalue",
                title = "Mov10_oe vs Irrel_kd",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 4.0)

# Perform variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Plot PCA
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = condition)) +
 geom_point(size = 3) +
 xlab(paste0("PC1: ", percentVar[1], "% variance")) +
 ylab(paste0("PC2: ", percentVar[2], "% variance")) +
 ggtitle("PCA of Mov10_oe vs Irrel_kd")



# ------------------------------------------------------------------------------
# Hierarchical Clustering Heatmap of 20 most differentially expressed genes

# Select the top 20 differentially expressed genes by p-value
top_genes <- head(order(res$padj), 20)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)

# Plot heatmap
pheatmap(mat, cluster_rows = TRUE, show_rownames = TRUE, 
         cluster_cols = TRUE, annotation_col = sample_info)


# Order by absolute log2 fold change and p-value
top_10_up <- head(res[order(res$log2FoldChange, -res$pvalue), ], 10)
top_10_down <- head(res[order(-res$log2FoldChange, res$pvalue), ], 10)

# Display top 10 upregulated genes
top_10_up

# Display top 10 downregulated genes
top_10_down


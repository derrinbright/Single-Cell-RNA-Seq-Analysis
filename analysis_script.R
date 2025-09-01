dir.create("D:/R_libs", recursive = TRUE, showWarnings = FALSE)
.libPaths("D:/R_libs")
.libPaths()  # The first path listed should now be D:/R_libs
install.packages(c("dplyr", "patchwork", "ggplot2", "RColorBrewer"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("SingleR", "celldex"))
install.packages("Seurat")

setwd("D:/single-cell-project")
file.rename("analysis_script.r.txt", "analysis_script.R")
setwd("D:\\single-cell-project\\raw_data")
getwd()

# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(celldex)
library(RColorBrewer)

# Set seed for reproducibility
set.seed(1234)

# Read 10X data from a directory
data_dir <- "D:/single-cell-project/raw_data"
raw_data <- Read10X(data.dir = data_dir)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = raw_data,
                                 project = "B16_Melanoma",
                                 min.cells = 3, # Genes detected in ≥3 cells
                                 min.features = 200) # Cells with ≥200 genes

# Add metadata (manually specify genotype)
seurat_obj$genotype <- "WT"
seurat_obj$sample <- "sample1"

# Basic summary
seurat_obj
# View metadata
head(seurat_obj@meta.data)

# Step 2: Quality Control 
# Calculate the percentage of mitochondrial genes per cell
# In mouse datasets, mitochondrial genes typically start with "mt-"
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Visualize key QC metrics per cell using violin plots:
# - nFeature_RNA: number of genes detected per cell
# - nCount_RNA: total RNA counts per cell
# - percent.mt: mitochondrial gene percentage per cell

Idents(seurat_obj) <- "orig.ident"

VlnPlot(seurat_obj,
        features = c("nFeature_RNA","nCount_RNA", "percent.mt"),
        pt.size = 0.1,
        ncol = 3)

table(Idents(seurat_obj))
unique(Idents(seurat_obj))  # Returns: "B16_Melanoma"



# Filter out low-quality cells based on:
# - Cells expressing fewer than 200 genes (likely dead)
# - Cells expressing more than 6000 genes (potential doublets)
# - Cells with >15% mitochondrial gene content (stressed/dead)
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 &
                       nFeature_RNA < 6000 &
                       percent.mt < 15)

# Step 3: Normalization and Feature Selection
# Normalize the data using log-normalization:
# Transforms raw UMI counts into normalized values to account for sequencing depth
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify 2000 most variable genes across all cells
# These genes will be used for downstream dimensionality reduction and clustering
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Visualize the variable features
# The top 10 variable genes are labeled for interpretability
var_plot <- VariableFeaturePlot(seurat_obj)
LabelPoints(plot = var_plot, points = head(VariableFeatures(seurat_obj), 20), repel = TRUE)


# Step 4: Dimensionality Reduction & Clustering
# Scale and center the expression values of each gene across all cells
# This step is essential for PCA and downstream linear modeling
seurat_obj <- ScaleData(seurat_obj)

# Run PCA to reduce dimensionality and identify major sources of variation
# npcs = number of principal components to compute (here we use 50)
seurat_obj <- RunPCA(seurat_obj, npcs = 50)

# Visualize the standard deviation of each PC to decide how many PCs to retain
# A sharp drop ("elbow") in the plot helps determine the optimal number of PCs
ElbowPlot(seurat_obj, ndims = 50)

# Based on the elbow plot, we select the first 15 principal components
pcs <- 15

# Identify the k-nearest neighbors for each cell based on the selected PCs
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pcs)
# Perform graph-based clustering of the cells
# Resolution affects the granularity: higher = more clusters
seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)
# Run UMAP for 2D visualization of the clusters
seurat_obj <- RunUMAP(seurat_obj, dims = 1:pcs)

# Plot UMAP clustering result
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("UMAP: Clustered Cells")


# Step 5: Cell Type Annotation 
# Loading a built-in reference dataset (Mouse RNA-seq) from celldex package
ref <- celldex::MouseRNAseqData()


# Automatically assign cell types using SingleR
# Compares each single cell to reference cell types using expression profile similarity
annotations <- SingleR(
  test = GetAssayData(seurat_obj, slot = "data"),
  ref = ref, 
  labels = ref$label.main
)

# Add predicted cell type labels to the Seurat object metadata
seurat_obj$celltype <- annotations$labels

# Plot UMAP with predicted cell types
DimPlot(seurat_obj, group.by = "celltype", label = TRUE, repel = TRUE) + 
  ggtitle("UMAP: Annotated Cell Types")

# Step:6 Find marker genes for all clusters
# Only return genes with higher expression in the cluster compared to others
markers <- FindAllMarkers(seurat_obj, 
                          only.pos = TRUE,           # keep only positive markers
                          min.pct = 0.25,            # expressed in at least 25% of cells in a cluster
                          logfc.threshold = 0.25)    # minimum log-fold change

# 🔹 Extract top 5 markers per cluster based on avg_log2FC
library(dplyr)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# 🔹 Save markers to CSV file
write.csv(top_markers, "cluster_markers.csv")

# 🔹 Perform differential expression between two specific clusters
# You can change 'ident.1' and 'ident.2' to other cluster IDs or cell types
de_genes <- FindMarkers(seurat_obj,
                        ident.1 = 0,
                        ident.2 = 1,
                        min.pct = 0.25,
                        logfc.threshold = 0.25)

# 🔹 View top DE genes
head(de_genes)

# 🔹 Save to CSV
write.csv(de_genes, "cluster0_vs_cluster1_DEGs.csv")

# Add a column for significance
de_genes$gene <- rownames(de_genes)
de_genes$significant <- ifelse(de_genes$p_val_adj < 0.05 & abs(de_genes$avg_log2FC) > 0.5, "Yes", "No")

# Plot using ggplot2
library(ggplot2)
ggplot(de_genes, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Cluster 0 vs Cluster 1",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")


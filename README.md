# Single-Cell Transcriptomic Profiling of the B16 Melanoma's Microenvironment

**Project Status: Completed**

---

## 1. Project Overview

This project is a complete computational analysis of a public single-cell RNA sequencing (scRNA-seq) dataset. The primary objective was to implement an end-to-end bioinformatics workflow in R to process raw gene expression data from a B16 melanoma tumor model, identify the distinct cell populations, and characterize their biological identities.

The project successfully transformed raw, high-dimensional count data into a clear map of the tumor's cellular landscape. The analysis identified and annotated **14 distinct cell types**, revealing that the tumor is not a uniform mass of cancer cells but a complex ecosystem of immune and stromal cells.

-   **Dataset:** [GSE110746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110746) from Gene Expression Omnibus (GEO)
-   **Model:** B16-F10 Murine (Mouse) Melanoma

---

## 2. Background and Rationale

### Why This Project?

A tumor is not just a collection of cancer cells. It's a complex **Tumor Microenvironment (TME)**, a bustling "neighborhood" of immune cells, structural cells (fibroblasts), and blood vessel cells. These neighboring cells determine whether a tumor grows, spreads, or responds to treatment.

Traditional sequencing (bulk RNA-seq) is like putting this entire neighborhood into a blender—you get an *average* gene expression profile, but you lose all the detail.

**Single-cell RNA-seq** is a revolutionary technique that solves this. It's like interviewing every "resident" (each cell) one-by-one. This allows us to build a detailed census of the tumor, identifying every cell type and understanding its unique function. The goal of this project was to perform this "cellular census" on a standard B16 melanoma model.

### About the Data

-   **B16 Melanoma Model:** This is a widely-used mouse model for studying aggressive skin cancer. It's the "crash test dummy" of melanoma research, allowing scientists to reliably study tumor growth and test new therapies.
-   **GSE110746 Dataset:** This is the raw count data from the original study, which used 10x Genomics technology to profile the cells within this tumor model.

---

## 3. Tools and Technologies

This entire analysis was conducted in **R**. The core of the workflow was built using a suite of powerful bioinformatics packages:

-   **Core Language:** R
-   **Primary Analysis Toolkit:** **Seurat** (for data handling, QC, normalization, clustering, and visualization)
-   **Cell Type Annotation:** **SingleR** & **celldex** (for automated cell identity prediction)
-   **Data Manipulation:** **dplyr** (for filtering and ranking genes)
-   **Data Visualization:** **ggplot2** & **patchwork** (for creating custom plots like the volcano plot)
-   **Data Source:** NCBI Gene Expression Omnibus (GEO)

---

## 4. Methodology: A Step-by-Step Workflow

The project was executed in a logical pipeline, moving from raw data to biological insight.

### Phase 1: Data Acquisition and Quality Control (QC)

The raw 10x Genomics data (`matrix.mtx`, `barcodes.tsv`, `features.tsv`) was loaded into R.

1.  **Create Seurat Object:** The data was loaded into a `Seurat Object`, the central container for the entire analysis.
2.  **Calculate QC Metrics:** We calculated the percentage of mitochondrial genes (`percent.mt`) for every cell. A high percentage can indicate a stressed or dying cell.
3.  **Visualize QC:** Before filtering, we created violin plots to visualize the distribution of three key metrics:
    -   `nFeature_RNA`: The number of unique genes detected in each cell.
    -   `nCount_RNA`: The total number of RNA molecules detected in each cell.
    -   `percent.mt`: The percentage of mitochondrial genes.
4.  **Filtering (Subsetting):** Based on the plots, we filtered out low-quality data. Cells with too few genes (likely empty droplets), too many genes (potential doublets), or high mitochondrial DNA were removed to ensure a clean dataset of healthy, single cells.

<img src="images/Violin_Plot.png" width="600"/>  

*Figure 1: QC violin plots showing the distribution of genes, RNA counts, and mitochondrial percentage before filtering.*

### Phase 2: Normalization and Feature Selection

1.  **Normalization (`NormalizeData`):** Raw UMI counts were normalized using a log-transformation. This step is crucial for comparing gene expression levels across cells that may have had different sequencing depths.
2.  **Find Variable Features (`FindVariableFeatures`):** We identified the 2,000 most highly variable genes across all cells. Focusing on these genes, which show the most biological variation, allows us to find the patterns that separate cell types.
3.  **Visualize Variable Features:** The variable feature plot shows the average expression versus the standardized variance. The red dots represent the 2,000 most variable genes (like `Sparc` and `Col1a2`) that will be used for the next steps.

<img src="images/Variable_Feature_Plot.png" width="500"/>

*Figure 2: A plot identifying the 2,000 most variable genes (red) that drive the biological differences in the dataset.*

### Phase 3: Dimensionality Reduction and Clustering

1.  **Scale Data (`ScaleData`):** The expression of the 2,000 variable genes was scaled to give each gene equal weight.
2.  **Run PCA (`RunPCA`):** We performed Principal Component Analysis (PCA) to reduce the 2,000-dimensional data into its main patterns of variation (we computed 50 PCs).
3.  **Select Dimensions:** An `ElbowPlot` (not shown) helped us determine that the first **15 PCs** captured the vast majority of the "signal" in the data.
4.  **Find Clusters (`FindNeighbors`, `FindClusters`):** Using these 15 PCs, we built a "social network" of cells, connecting cells that were most similar. The algorithm then identified "communities" or "cliques" within this network. These are the cell clusters.
5.  **Run UMAP (`RunUMAP`):** To visualize these clusters, we used UMAP to create a simple 2D "map" of the cells, where each "island" represents a distinct cluster.

---

## 5. Results and Analysis

This phase focused on interpreting the clusters and understanding their biological meaning.

### 5.1 Cell Type Annotation

The `FindClusters` function gave us abstract labels (Cluster 0, Cluster 1, etc.). To identify them, we used **SingleR**.

1.  **Annotation:** SingleR compares the gene expression profile of each cluster to a reference database of known mouse cell types (the `MouseRNAseqData` from celldex).
2.  **Labeling:** It then assigns the most likely cell type label to each cluster.
3.  **Visualization:** We then plotted the UMAP again, but this time colored by the new biological identities.

<img src="images/Annotated_Cell_Types.png" width="700"/>

*Figure 3: The final UMAP plot, where each cluster is annotated with its predicted cell type. This map clearly shows 14 distinct cell populations, including T cells, Macrophages, Monocytes, and Fibroblasts.*

This plot is the central finding of the project. It clearly demonstrates the cellular heterogeneity of the tumor microenvironment.

### 5.2 Differential Expression (DE) Analysis

To understand *why* these cell types are different, we identified their unique gene signatures.

1.  **Find All Markers (`FindAllMarkers`):** We performed DE analysis for *each* cluster to find the genes that define it. The top 5 markers for each cluster were saved to a CSV file.
2.  **Compare Two Clusters (`FindMarkers`):** We also performed a direct comparison between two specific clusters (Cluster 0 and Cluster 1) to find all genes that were differentially expressed between them.
3.  **Visualize DEGs:** We visualized these results using a **Volcano Plot**.

<img src="images/Volcano_Plot.png" width="600"/>

*Figure 4: A volcano plot comparing Cluster 0 vs. Cluster 1. The X-axis is the log-fold change (how much more/less expressed a gene is). The Y-axis is the statistical significance. The red dots represent genes that are both highly significant and show a large change, making them key biological differentiators between these two cell populations.*

---

## 6. Conclusion

This project successfully implemented a complete scRNA-seq analysis workflow in R. Starting from raw count data, the pipeline was able to perform quality control, normalize the data, and use unsupervised clustering to identify the underlying cell populations.

The key takeaway is that the B16 melanoma model is not a uniform mass but a complex ecosystem. The analysis revealed a rich cellular landscape of **14 distinct cell types**, including a diverse infiltrate of immune and stromal cells. This annotated "cell atlas" provides a foundational map of the tumor microenvironment, which is the essential first step for any future study into how these different cell types communicate and how the tumor might respond to therapy.

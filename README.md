## Overview

Understanding cell-cell interactions is pivotal for unraveling the complexities of biological processes such as development, immune responses, and cancer progression. The temporal dynamics of these interactions can offer significant insights into how cells communicate and influence each other over time. To facilitate this analysis, two specialized functions have been developed, each addressing a critical aspect of temporal cell-cell communication:

## `TemporalCCI1`

This function is designed to calculate the temporal correlations between ligands expressed by one cell type and receptors expressed by another cell type. By analyzing how these interactions evolve across a defined pseudotime trajectory, `TemporalCCI1` helps researchers identify dynamic patterns in cell communication. This approach is particularly useful in studies that require an understanding of the timing and evolution of cell-cell interactions, such as developmental biology or cancer progression research.

## `TemporalCCI2`

`TemporalCCI2` extends the analysis by inferring the impact of one cell type on the evolutionary process of another based on their communication. The function divides the pseudotime into discrete intervals and utilizes the CellChat framework to analyze ligand-receptor interactions across these stages. By doing so, `TemporalCCI2` provides insights into how one cell type might influence the development or differentiation of another, helping to identify critical communication pathways that drive cellular evolution.

Together, these functions offer a comprehensive toolkit for exploring the temporal dynamics of cell-cell interactions, enabling a more nuanced understanding of how cellular communication shapes biological processes over time.


## Installation
```r
# install.packages("devtools")
devtools::install_github('Zaoqu-Liu/TimeCCI')
```

## TemporalCCI1 Tutorial
This function performs a temporal cell-cell interaction (CCI) analysis between two specified cell types over pseudotime.
It leverages ligand-receptor interaction data to assess the correlation between ligand expression in one cell type
and receptor expression in another across a defined pseudotime trajectory.

```r
# Clear workspace and load necessary libraries
rm(list = ls())
library(Seurat)
library(tidyverse)
library(monocle)

# Load the Seurat object
sc <- readRDS('example.rds')

# Preprocess the data: normalization, feature selection, scaling, PCA, and UMAP
# - NormalizeData: Normalizes the expression data for each gene.
# - FindVariableFeatures: Identifies the most variable genes for downstream analysis.
# - ScaleData: Scales the data for PCA.
# - RunPCA: Performs Principal Component Analysis for dimensionality reduction.
# - RunUMAP: Further reduces the data to 2D space for visualization.
sc <- sc %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 500) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(reduction = 'pca', dims = 1:30)

# Set cell identities based on the cell type annotation in metadata
Idents(sc) <- sc$Cell.Type

# Visualize the UMAP plot to see the clustering of cells
UMAPPlot(sc)

# Convert Seurat data to Monocle's format for pseudotime analysis
expr_matrix <- as(as.matrix(sc@assays$RNA@counts), 'sparseMatrix')
dim(expr_matrix)

# Extract metadata and prepare feature data for Monocle
p_data <- sc@meta.data 
p_data$celltype <- sc@active.ident
f_data <- data.frame(gene_short_name = row.names(sc), row.names = row.names(sc))

# Create AnnotatedDataFrame objects for phenodata and featuredata
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)

# Create the CellDataSet object for Monocle
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

# Estimate size factors and dispersions for normalization
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# Perform differential gene expression testing
diff_test_res <- differentialGeneTest(
  cds,
  fullModelFormulaStr = "~celltype",
  cores = 5
)

# Filter for significant genes with q-value < 0.01
diff_test_res1 <- diff_test_res[diff_test_res$qval < 0.01,]

# Select top 500 genes for ordering cells
ordergenes <- row.names(diff_test_res1)[1:500]
cds2 <- setOrderingFilter(cds, ordering_genes = ordergenes)

# Reduce dimensionality and order cells using DDRTree
cds2 <- reduceDimension(cds2, method = 'DDRTree', max_components = 2, num_dim = 10, cores = 13)
cds2 <- orderCells(cds2)

# Plot the cell trajectory colored by pseudotime and cell type
plot_cell_trajectory(cds2)
plot_cell_trajectory(cds2, color_by = 'celltype')

# Add pseudotime data to the Seurat object
sc$Pseudotime <- cds2$Pseudotime

# Perform temporal cell-cell interaction analysis using TemporalCCI1
tmp1 <- TemporalCCI1(sc)
```

## TemporalCCI2 Tutorial
This function first divides the pseudotime into discrete intervals and assigns each cell to an interval. It then uses
the CellChat package to analyze cell-cell communication based on ligand-receptor interactions, computing the communication
probabilities for each interaction across the pseudotime intervals. The function also performs correlation analyses
between these probabilities and the pseudotime stages, providing insights into how interactions evolve over time.

```r
# Clear the workspace and start a new analysis
rm(list = ls())
library(Seurat)
library(tidyverse)
library(monocle)

# Load another Seurat object for analysis
sc <- seu.example

# Subset the data into two groups: epithelial cells (EpiT) and mCAF cells
sc_epit <- subset(sc, Cell.Type == 'EpiT')
sc_mcaf <- subset(sc, Cell.Type == 'mCAF')

# Preprocess the mCAF subset: normalization, feature selection, scaling, PCA, and UMAP
sc_mcaf <- sc_mcaf %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 500) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(reduction = 'pca', dims = 1:30)

# Convert the mCAF data to Monocle's format for pseudotime analysis
expr_matrix <- as(as.matrix(sc_mcaf@assays$RNA@counts), 'sparseMatrix')
dim(expr_matrix)

# Extract metadata and prepare feature data for Monocle
p_data <- sc_mcaf@meta.data 
p_data$celltype <- sc_mcaf@active.ident
f_data <- data.frame(gene_short_name = row.names(sc_mcaf), row.names = row.names(sc_mcaf))

# Create AnnotatedDataFrame objects for phenodata and featuredata
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)

# Create the CellDataSet object for Monocle
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

# Estimate size factors and dispersions for normalization
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# Identify highly variable genes
disp_table <- dispersionTable(cds)
disp.gene <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= dispersion_fit)$gene_id

# Perform dimensionality reduction and order cells using these genes
cds2 <- setOrderingFilter(cds, ordering_genes = disp.gene)
cds2 <- reduceDimension(cds2, method = 'DDRTree', max_components = 2, num_dim = 10, cores = 13)
cds2 <- orderCells(cds2)

# Plot the pseudotime trajectory for mCAF cells
plot_cell_trajectory(cds2)

# Add pseudotime data to the mCAF Seurat object
sc_mcaf$Pseudotime <- cds2$Pseudotime

# Merge the mCAF and epithelial cell subsets back together
sc <- merge(sc_epit, sc_mcaf)

# Preprocess the merged dataset: normalization, feature selection, scaling, PCA, and UMAP
sc <- sc %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 500) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(reduction = 'pca', dims = 1:30)

# Set cell identities based on the cell type annotation in metadata
Idents(sc) <- sc$Cell.Type

# Visualize the UMAP plot for the merged dataset
UMAPPlot(sc)

# Perform temporal cell-cell interaction analysis using TemporalCCI2
tmp2 <- TemporalCCI2(sc)
```

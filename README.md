



# TemporalCCI Usage Tutorial

This tutorial demonstrates how to use the `TemporalCCI1` function with an example Seurat object. We will preprocess the data, perform dimensionality reduction, calculate pseudotime using Monocle, and analyze cell-cell interactions over pseudotime using `TemporalCCI1`.

## Complete Code

```r
# Clear workspace and load necessary libraries
rm(list = ls())
library(Seurat)
library(tidyverse)
library(monocle)

# Load the Seurat object (example)
sc <- seu.example

# Preprocessing steps: normalization, feature selection, scaling, PCA, and UMAP
sc <- sc %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 500) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(reduction = 'pca', dims = 1:30)

# Set cell identities
Idents(sc) <- sc$Cell.Type

# Visualize the UMAP plot
UMAPPlot(sc)

# Convert Seurat data to Monocle's format
expr_matrix <- as(as.matrix(sc@assays$RNA@counts), 'sparseMatrix')
dim(expr_matrix)

# Extract metadata and prepare feature data
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

# Perform temporal cell-cell interaction analysis
tmp1 <- TemporalCCI1(sc)

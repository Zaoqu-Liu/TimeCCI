#' Temporal Cell-Cell Interaction Analysis via Cell Alignment in Pseudotime
#'
#' This function performs a temporal cell-cell interaction (CCI) analysis between two specified cell types over pseudotime.
#' It leverages ligand-receptor interaction data to assess the correlation between ligand expression in one cell type
#' and receptor expression in another across a defined pseudotime trajectory.
#'
#' @param seu.obj A Seurat object containing single-cell RNA-seq data and metadata.
#' @param CellType.name A string specifying the name of the cell type annotation in the Seurat object metadata. Default is "Cell.Type".
#' @param Pseudotime.name A string specifying the name of the pseudotime annotation in the Seurat object metadata. Default is "Pseudotime".
#' @param CCI.database A data frame containing ligand-receptor interaction pairs with columns "interaction_name", "ligand.symbol", and "receptor.symbol".
#' @param source.celltype A string specifying the source cell type for the ligand in the CCI analysis. Default is "mCAF".
#' @param target.celltype A string specifying the target cell type for the receptor in the CCI analysis. Default is "EpiT".
#' @param window.size A numeric value defining the window size for trajectory alignment. Default is 0.1.
#' @param num.pts An integer specifying the number of points for trajectory interpolation. Default is 200.
#' @param num.cores An integer specifying the number of CPU cores to use for parallel processing. Default is 10.
#' @param granger.cutoff A numeric value for the p-value cutoff in Granger causality tests. Default is 0.01.
#'
#' @return A list containing:
#' \item{res}{A data frame with correlation results between ligand and receptor pairs, including Pearson and Spearman correlation coefficients and p-values.}
#' \item{source.align.matrix}{A data frame with the aligned and interpolated expression matrix of the source cell type.}
#' \item{target.align.matrix}{A data frame with the aligned and interpolated expression matrix of the target cell type.}
#'
#' @details
#' This function is designed to analyze temporal changes in cell-cell communication by evaluating ligand-receptor pair correlations
#' over pseudotime. It uses interpolation methods to align gene expression trajectories, enabling the comparison of expression patterns
#' between a source and a target cell type. The function filters out ligand-receptor pairs that are not present in the RNA assay of
#' the Seurat object. Correlation analyses are performed using both Pearson and Spearman methods, and the results are compiled into
#' a comprehensive data frame. Parallel processing is supported to expedite computations.
#'
#' @export
#' @author
#' Zaoqu Liu; Email: liuzaoqu@163.com
TemporalCCI1 <- function(
    seu.obj,
    CellType.name = "Cell.Type",
    Pseudotime.name = "Pseudotime",
    CCI.database = LRdf,
    source.celltype = "mCAF",
    target.celltype = "EpiT",
    window.size = 0.1,
    num.pts = 200,
    num.cores = 10,
    granger.cutoff = 0.01) {
  # Remove 'cell_id' column from metadata to prevent duplication
  seu.obj@meta.data <- seu.obj@meta.data[, colnames(seu.obj@meta.data) != "cell_id"]

  # Define a function for Min-Max scaling
  MinMaxScale <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }

  # Define a function to replace NA values with 0 in a data frame
  myRemoveNA <- function(df) {
    df[is.na(df)] <- 0
    return(df)
  }

  # Extract ligand-receptor pairs and their respective gene lists from the database
  LRpairs <- CCI.database[, 1]
  Lgenelist <- CCI.database[, 2]
  Rgenelist <- CCI.database[, 3]

  # Use data from the RNA assay
  tmp.data <- Seurat::GetAssayData(object = seu.obj, slot = "data", assay = "RNA")
  gene_symbols <- rownames(tmp.data)

  # Filter out ligand-receptor pairs that are not present in the RNA data
  index.remove <- which(!(Lgenelist %in% gene_symbols) | !(Rgenelist %in% gene_symbols))
  LRpairs <- LRpairs[-index.remove]
  Lgenelist <- Lgenelist[-index.remove]
  Rgenelist <- Rgenelist[-index.remove]

  tmp.df <- seu.obj[[]]

  # Check if the specified cell type and pseudotime columns exist in the metadata
  if (!CellType.name %in% colnames(tmp.df)) {
    stop("Please add CellType annotation in Seurat object!")
  }

  if (!Pseudotime.name %in% colnames(tmp.df)) {
    stop("Please add Pseudotime information in seurat object!")
  }

  # Rename columns for easier reference
  colnames(tmp.df)[colnames(tmp.df) == CellType.name] <- "celltype"
  colnames(tmp.df)[colnames(tmp.df) == Pseudotime.name] <- "pseudotime"

  # Extract and scale metadata for the source and target cell types
  tmp.cell.meta.1 <- tmp.df %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::filter(celltype == source.celltype) %>%
    dplyr::arrange(pseudotime)
  tmp.cell.meta.1$pseudotime <- MinMaxScale(tmp.cell.meta.1$pseudotime)

  tmp.cell.meta.2 <- tmp.df %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::filter(celltype == target.celltype) %>%
    dplyr::arrange(pseudotime)
  tmp.cell.meta.2$pseudotime <- MinMaxScale(tmp.cell.meta.2$pseudotime)

  # Extract expression matrices for the source and target cell types
  tmp.mat.1 <- tmp.data[, tmp.cell.meta.1$cell_id]
  tmp.mat.2 <- tmp.data[, tmp.cell.meta.2$cell_id]

  tmp.mat.pseudotime.1 <- tmp.cell.meta.1$pseudotime
  tmp.mat.pseudotime.2 <- tmp.cell.meta.2$pseudotime

  names(tmp.mat.pseudotime.1) <- tmp.cell.meta.1$cell_id
  names(tmp.mat.pseudotime.2) <- tmp.cell.meta.2$cell_id

  # Interpolate and scale the expression data along the pseudotime trajectory
  inter.tmp.mat.1 <- cellAlign::interWeights(
    expDataBatch = tmp.mat.1,
    trajCond = tmp.mat.pseudotime.1,
    winSz = window.size,
    numPts = num.pts
  )
  inter.tmp.mat.2 <- cellAlign::interWeights(
    expDataBatch = tmp.mat.2,
    trajCond = tmp.mat.pseudotime.2,
    winSz = window.size,
    numPts = num.pts
  )

  inter.tmp.mat.1 <- cellAlign::scaleInterpolate(inter.tmp.mat.1)
  inter.tmp.mat.2 <- cellAlign::scaleInterpolate(inter.tmp.mat.2)
  time <- inter.tmp.mat.1$traj

  inter.tmp.mat.1 <- myRemoveNA(inter.tmp.mat.1$scaledData)
  inter.tmp.mat.2 <- myRemoveNA(inter.tmp.mat.2$scaledData)

  # Use parallel processing to calculate correlations for each ligand-receptor pair
  future::plan("multisession", workers = num.cores)
  res <- furrr::future_map_dfr(seq_along(Lgenelist), \(ii){
    x <- inter.tmp.mat.1[Lgenelist[ii], ]
    y <- inter.tmp.mat.2[Rgenelist[ii], ]
    pp <- stats::cor.test(x, y, method = "pearson")
    ss <- stats::cor.test(x, y, method = "spearman")
    return(
      data.frame(
        LRpair = LRpairs[ii], Ligand = Lgenelist[ii], Receptor = Rgenelist[ii],
        PCC = ifelse(is.na(pp$estimate), 0, pp$estimate),
        PCP = ifelse(is.na(pp$p.value), 1, pp$p.value),
        SCC = ifelse(is.na(ss$estimate), 0, ss$estimate),
        SCP = ifelse(is.na(ss$p.value), 1, ss$p.value),
        row.names = FALSE
      )
    )
  }, .progress = TRUE)

  return(list(
    res = res,
    source.align.matrix = as.data.frame(t(inter.tmp.mat.1)),
    target.align.matrix = as.data.frame(t(inter.tmp.mat.2))
  ))
}

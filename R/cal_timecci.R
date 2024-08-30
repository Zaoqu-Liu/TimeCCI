#' Temporal Cell-Cell Interaction Analysis in Single-Cell RNA-seq Data
#'
#' This function performs temporal cell-cell interaction analysis in single-cell RNA-seq data using the CellChat framework.
#' It calculates correlations between ligand-receptor interaction probabilities and pseudotime across different cell states.
#'
#' @param seu.obj A Seurat object that contains single-cell RNA-seq data. The object must include metadata columns for pseudotime and cell type.
#' @param Pseudotime.name The name of the metadata column in `seu.obj` that contains pseudotime information. Default is 'Pseudotime'.
#' @param CellType.name The name of the metadata column in `seu.obj` that contains cell type information. Default is 'CellType'.
#' @param discrete.numer The number of discrete time points to divide the pseudotime into. Default is 20.
#' @param min.cells The minimum number of cells required for a cell-cell interaction to be considered. Default is 10.
#' @param CCI.database The cell-cell interaction database to use. Default is `CellChat::CellChatDB.human`.
#' @param interaction.name The name of the column in the CellChat LRsig table that contains the interaction name. Default is "interaction_name_2".
#' @return A list containing:
#' \item{cor.res}{A data frame with the correlation results, including p-values and correlation coefficients.}
#' \item{cal.res}{A data frame with the calculated ligand-receptor interaction probabilities over time.}
#' \item{cellchat}{The updated CellChat object after the analysis.}
#' @details
#' The function first removes any missing data from the Seurat object's metadata. It then divides the pseudotime into discrete intervals
#' and adds this information as new metadata to the Seurat object. The CellChat object is then created and used to compute communication
#' probabilities between cells based on their temporal states. Finally, correlations between ligand-receptor interaction probabilities
#' and pseudotime are calculated.
#' @export
cal_timecci <- function(seu.obj,
                        Pseudotime.name = "Pseudotime",
                        CellType.name = "CellType",
                        discrete.numer = 20,
                        min.cells = 10,
                        CCI.database = CellChat::CellChatDB.human,
                        interaction.name = "interaction_name_2") {
  # Remove missing data from the Seurat object's metadata
  df <- stats::na.omit(seu.obj@meta.data)

  # Divide pseudotime into discrete intervals and label them
  df$CellTime <- cut(df[[Pseudotime.name]],
    breaks = quantile(df[[Pseudotime.name]], probs = seq(0, 1, by = 1 / discrete.numer)),
    include.lowest = TRUE,
    labels = paste0("T", 1:discrete.numer)
  )
  df$CellTime <- paste0("Z-", as.character(df$CellTime))

  # Add the modified metadata back to the Seurat object
  seu.obj <- Seurat::AddMetaData(seu.obj, df)
  seu.obj$CellTime2 <- ifelse(is.na(seu.obj$CellTime), "A", seu.obj$CellTime)
  Seurat::Idents(seu.obj) <- seu.obj$CellTime2

  # Extract RNA counts and metadata
  x <- seu.obj@assays$RNA@counts
  data <- list(data = x, meta = seu.obj@meta.data)
  data.input <- data$data
  meta <- data$meta

  # Create CellChat object and assign the appropriate database
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "CellTime2")
  cellchat@DB <- CCI.database

  # Perform CellChat analysis steps
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, trim = 0.1, nboot = 10, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)

  # Extract communication probability data
  df <- cellchat@net[["prob"]]
  LR <- cellchat@LR[["LRsig"]]
  LR <- LR[match(dimnames(df)[[3]], LR$interaction_name), ]

  # Prepare dataframe for correlation analysis
  column_names <- rownames(df[, , 1])[-1]
  LR[column_names] <- NA
  for (i in LR$interaction_name) {
    LR[LR$interaction_name == i, (ncol(cellchat@LR[["LRsig"]]) + 1):ncol(LR)] <- df[1, -1, i]
  }

  # Transform data to long format for correlation analysis
  LR2 <- LR[, c(which(colnames(LR) == interaction.name), (ncol(cellchat@LR[["LRsig"]]) + 1):ncol(LR))]
  LR2 <- tidyr::pivot_longer(LR2, cols = 2:ncol(LR2), names_to = "CellTime", values_to = "Prob")
  LR2$CellTime2 <- as.numeric(gsub("Z-T", "", LR2$CellTime))
  LR2 <- as.data.frame(LR2)

  # Correlation analysis between interaction probability and pseudotime
  res <- data.frame()
  for (i in unique(LR2[, interaction.name])) {
    dd <- LR2[LR2[, interaction.name] == i, ]
    fit <- cor.test(dd$Prob, dd$CellTime2)
    res <- rbind(res, data.frame(name = i, pval = fit$p.value, scc = fit$estimate))
  }

  rownames(res) <- NULL

  # Return results
  return(list(
    cor.res = res,
    cal.res = LR2,
    cellchat = cellchat
  ))
}

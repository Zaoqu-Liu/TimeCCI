#' Temporal Cell-Cell Interaction Analysis with Discrete Pseudotime Binning
#'
#' This function performs cell-cell interaction (CCI) analysis over discrete pseudotime intervals.
#' It utilizes the CellChat framework to identify and analyze ligand-receptor interactions across different
#' stages of pseudotime, and computes correlations between interaction probabilities and pseudotime stages.
#'
#' @param seu.obj A Seurat object containing single-cell RNA-seq data and metadata.
#' @param Pseudotime.name A string specifying the name of the pseudotime annotation in the Seurat object metadata. Default is "Pseudotime".
#' @param CellType.name A string specifying the name of the cell type annotation in the Seurat object metadata. Default is "Cell.Type".
#' @param discrete.numer An integer specifying the number of discrete pseudotime intervals to create. Default is 20.
#' @param CCI.database A CellChat database object containing ligand-receptor interactions. Default is `CellChat::CellChatDB.human`.
#' @param interaction.name A string specifying the interaction name column in the CCI database. Default is "interaction_name_2".
#'
#' @return A list containing:
#' \item{cor.res}{A data frame with correlation results between interaction probabilities and pseudotime stages, including Pearson and Spearman correlation coefficients and p-values.}
#' \item{cal.res.wide}{A wide-format data frame with interaction probabilities across pseudotime stages.}
#' \item{cal.res.long}{A long-format data frame with interaction probabilities and corresponding pseudotime stages.}
#' \item{cellchat}{The CellChat object used for the analysis, containing all intermediate results.}
#'
#' @details
#' This function first divides the pseudotime into discrete intervals and assigns each cell to an interval. It then uses
#' the CellChat package to analyze cell-cell communication based on ligand-receptor interactions, computing the communication
#' probabilities for each interaction across the pseudotime intervals. The function also performs correlation analyses
#' between these probabilities and the pseudotime stages, providing insights into how interactions evolve over time.
#'
#' @export
#' @author
#' Zaoqu Liu; Email: liuzaoqu@163.com
TemporalCCI2 <- function(seu.obj,
                         Pseudotime.name = "Pseudotime",
                         CellType.name = "Cell.Type",
                         discrete.numer = 20,
                         CCI.database = CellChat::CellChatDB.human,
                         interaction.name = "interaction_name_2") {
  # Remove missing data from the Seurat object's metadata
  df <- stats::na.omit(seu.obj@meta.data)

  # Divide pseudotime into discrete intervals and label them
  df$CellTime <- cut(df[[Pseudotime.name]],
    breaks = stats::quantile(df[[Pseudotime.name]], probs = seq(0, 1, by = 1 / discrete.numer)),
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
  cellchat <- CellChat::createCellChat(object = data.input, meta = meta, group.by = "CellTime2")
  cellchat@DB <- CCI.database

  # Perform CellChat analysis steps
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, raw.use = TRUE, trim = 0.1, nboot = 10, population.size = TRUE)

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
    pp <- stats::cor.test(dd$Prob, dd$CellTime2, method = "pearson")
    ss <- stats::cor.test(dd$Prob, dd$CellTime2, method = "spearman")
    res <- rbind(res, data.frame(
      LRpair = i,
      PCC = ifelse(is.na(pp$estimate), 0, pp$estimate),
      PCP = ifelse(is.na(pp$p.value), 1, pp$p.value),
      SCC = ifelse(is.na(ss$estimate), 0, ss$estimate),
      SCP = ifelse(is.na(ss$p.value), 1, ss$p.value),
      row.names = FALSE
    ))
  }

  # Return results
  return(list(
    cor.res = res,
    cal.res.wide = LR,
    cal.res.long = LR2,
    cellchat = cellchat
  ))
}

#' LRdf: Ligand-Receptor Interaction Database
#'
#' A data frame containing ligand-receptor interactions derived from the CellChat database.
#' This dataset includes the interaction name, ligand symbol, and receptor symbol.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{interaction_name}{The name of the ligand-receptor interaction.}
#'   \item{ligand.symbol}{The symbol of the ligand gene.}
#'   \item{receptor.symbol}{The symbol of the receptor gene.}
#' }
#' @source CellChatDB.human$interaction
#'
#' @examples
#' data(LRdf)
#' head(LRdf)
"LRdf"

#' seu.example: Seurat Object for mCAF and Tumor Epithelial Cells
#'
#' A Seurat object containing single-cell RNA-seq data and metadata for mCAF and tumor epithelial cells.
#' This object includes the RNA assay and relevant metadata such as cell type and pseudotime annotations.
#'
#' @format A Seurat object with RNA assay data and metadata:
#' \describe{
#'   \item{RNA}{The RNA assay data containing gene expression values.}
#'   \item{meta.data}{Metadata including cell type annotations and pseudotime values.}
#' }
#'
#' @examples
#' data(seu.example)
#' seu.example
"seu.example"

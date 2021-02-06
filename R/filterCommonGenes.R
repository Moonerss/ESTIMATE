#' @name filterCommonGenes
#' @title Intersect input data with 10,412 common genes
#' @description This function unifies different number of genes per platform against 10,412 common genes.
#' @param expr a matrix including gene-level expression data, rownames is genes, colnames is samples.
#' @param id character string indicating which gene identifier to use when matching.
#' @details The number of genes in expression data is different for each platform and this difference influences the computational results of stromal and immune scores. To compare stromal, immune and
#'   ESTIMATE scores across different platforms or calculate ESTIMATE-based tumor purity using
#'   Affymetrix expression data, users need to unify the gene identifiers of the input data against the
#'   common genes.
#' @return return a \emph{data.frame} of the result.
#' @author Erjie Zhao <2055469819@qq.com>
#' @export
#' @examples
#'   file <- system.file("extdata", "sample_input.txt", package="ESTIMATE")
#'   expression <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
#'   merged_expr <- filterCommonGenes(expr = data.matrix(expression), id = "GeneSymbol")
filterCommonGenes <- function(expr, id=c("GeneSymbol", "EntrezID")) {
  ## Check arguments
  stopifnot(is.matrix(expr) && is.numeric(expr))
  id <- match.arg(id)

  expr <- as.data.frame(expr)
  merged.df <- merge(common_genes, expr, by.x=id, by.y="row.names")
  rownames(merged.df) <- merged.df$GeneSymbol
  merged.df <- merged.df[, -1:-ncol(common_genes)]
  print(sprintf("Merged dataset includes %d genes (%d mismatched).",
                nrow(merged.df),
                nrow(common_genes) - nrow(merged.df)))

  ## get the object
  res <- data.frame(NAME = rownames(merged.df), Description=rownames(merged.df), merged.df)
  rownames(res) <- NULL
  res <- res[order(res$NAME),]
  return(res)
}

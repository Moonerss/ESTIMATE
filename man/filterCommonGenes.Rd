% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterCommonGenes.R
\name{filterCommonGenes}
\alias{filterCommonGenes}
\title{Intersect input data with 10,412 common genes}
\usage{
filterCommonGenes(expr, id = c("GeneSymbol", "EntrezID"))
}
\arguments{
\item{expr}{a matrix including gene-level expression data, rownames is genes, colnames is samples.}

\item{id}{character string indicating which gene identifier to use when matching.}
}
\value{
return a \emph{data.frame} of the result.
}
\description{
This function unifies different number of genes per platform against 10,412 common genes.
}
\details{
The number of genes in expression data is different for each platform and this difference influences the computational results of stromal and immune scores. To compare stromal, immune and
  ESTIMATE scores across different platforms or calculate ESTIMATE-based tumor purity using
  Affymetrix expression data, users need to unify the gene identifiers of the input data against the
  common genes.
}
\examples{
  file <- system.file("extdata", "sample_input.txt", package="ESTIMATE")
  expression <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  merged_expr <- filterCommonGenes(expr = data.matrix(expression), id = "GeneSymbol")
}
\author{
Erjie Zhao <2055469819@qq.com>
}

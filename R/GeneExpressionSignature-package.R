#' Introduction to the GeneExpressionSignature Package
#' 
#' 
#' The \pkg{GeneExpressionSignature} add-on is an implementation of computing
#' distances among preprocessed gene-expression profiles of samples for .  The
#' distances can be used to detect similarities among the signatures of drugs,
#' diseases, and biological states of interest, and construct connectivity map.
#' 
#' This package contains functions for the distances computation based on gene
#' expression signature. First, list of genes is ranked according to their
#' expression ratios to produce the Prototype Ranked List (PRL). Second, all
#' the PRLs with the same state are aggregated by \code{\link{RankMerging}}
#' functions. Finally, all the ranked lists are made as one input of the
#' \code{\link{ScorePGSEA}} and \code{\link{ScoreGSEA}} functions to compute
#' the pairwise distances.
#' @name GeneExpressionSignature_package
NULL


#' sample data, a subset of the C-MAP
#' 
#' sample data, a subset of the C-MAP as , which is a collection of 50
#' genome-wide transcriptional expression data from cultured human cells
#' treated with 15 different small molecules.
#' 
#' 
#' @name exampleSet
#' @docType data
#' @format A ExpressionSet: assay data represents the 50 genome-wide
#' transcriptional expression data, phenotypic data describes 15 different
#' small molecules corresponds to the expression data (assay data).
#' @references \url{http://www.sciencemag.org/content/313/5795/1929.short} Lamb
#' et al., The Connectivity Map: Using Gene-Expression Signatures to Connect
#' Small Molecules, Genes, and Disease, science 2006
#' @keywords data datasets
#' @examples
#' 
#' data(exampleSet)
#' 
NULL









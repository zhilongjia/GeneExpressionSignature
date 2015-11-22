#' ScoreX
#' 
#' 
#' ScoreX: The eXtreme simlarity score
#' 
#' @param obj Either a numeric (logFC) matrix, a data.frame, or an ExpressionSet object.
#' @param metric simlarity metric. Such as XSum, XCos.
#' @param SignatureLength the length of "gene signature".
#' @param scale scale or not.
#' @param ... other parameters
#' @import Biobase
#' @export
#' @author Zhilong Jia
#' @references 
#' Cheng, Jie, and Lun Yang. "Comparing gene expression similarity metrics for 
#' connectivity map." Bioinformatics and Biomedicine (BIBM), 2013 IEEE 
#' International Conference on. IEEE, 2013.
#' 
#' @examples 
#' data(exampleSet)
#' ScoreXSum <- ScoreX(logfcData)
#' ScoreXCos <- ScoreX(logfcData, metric="XCos")

ScoreX <- function(obj, metric="XSum", SignatureLength=500, scale=TRUE, ...) {
    
    # Check parameter
    switch(class(obj),
           matrix = mat <- as.data.frame(obj),
           ExpressionSet = mat <- Biobase::exprs(obj),
           data.frame = { mat <- obj },
           stop("argument 'obj' must be a matrix, data.frame, or ExpressionSet 
            object")
    )
    SignatureLength <- round(SignatureLength)
    if (SignatureLength <1 | SignatureLength> nrow(mat)) {
        stop("SignatureLength should be more than 0 and less than the number of genes.")
    }
    
    # Top genes for each sample
    geneUp <- apply(mat, 2, function(x){rownames(mat)[order(x, decreasing = TRUE)[1:SignatureLength]] })
    geneDown <- apply(mat, 2, function(x){rownames(mat)[order(x, decreasing = FALSE)[1:SignatureLength]] })
    
    # Scoring
    score_mat <- matrix(0, ncol=ncol(mat), nrow=ncol(mat), dimnames=list(colnames(mat) , colnames(mat)))
    
    if (metric == "XSum") {
        for (i in colnames(mat)){
            for (j in colnames(mat)) {
                
                overLap_up <- intersect(geneUp[,i], c(geneUp[,j], geneDown[,j]) )
                overLap_down <- intersect(geneDown[,i], c(geneUp[,j], geneDown[,j]) )
                if (length(c(overLap_up, overLap_down)) != 0) {
                    score_xsum <- sum( mat[ overLap_up, j], na.rm=TRUE ) - sum( mat[ overLap_down, j], na.rm=TRUE)
                } else {
                    score_xsum <- 0
                }
                score_mat[i, j] <- score_xsum
            }
        }
    }
    
    if (metric=="XCos") {
        for (i in colnames(mat)){
            for (j in colnames(mat)) {
                overLap_gene <- intersect(union(geneUp[,i], geneDown[,i]), union(geneUp[,j], geneDown[,j]) )
                if (length(overLap_gene) >= 0) {
                    score_xcos <- crossprod(mat[overLap_gene, i], mat[overLap_gene, j] )
                } else {
                    score_xcos <- 0
                }
                score_mat[i, j] <- score_xcos
            }
        } 
    }
    # Normalization
    if (isTRUE(scale)) {
        score_mat <- apply(score_mat, 2, function(x){x/max(abs(x))})
    }
    
    return (score_mat)

}


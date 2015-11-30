#' Convert data.frame into ExpressionSet Object.
#' 
#' Convert data.frame into ExpressionSet Object.
#' @param dat a data.frame with gene in row, sample in column
#' @param state sample state. For multi instances.
#' @return an ExpressionSet object
#' @author Zhilong Jia
#' @import Biobase
#' @export
#' @examples 
#' dat <- data.frame(d1=1:600, D1=600:1, d2=sample(600))
#' datExp <- makeExpressionSet(dat, state=c("d", "D", "d"))
#' 
#' # Only for rankMerging
#' MergingSet <- RankMerging(datExp,"Spearman",weighted=TRUE)
#' #ds <- ScoreGSEA(datExp, 250)
#' #ds1 <- ScorePGSEA(datExp, 250)
#' #ds2 <- ScoreGSEA(MergingSet, 250)
#' #ds3 <- ScorePGSEA(MergingSet, 250)
#' 

makeExpressionSet <- function(dat, state=colnames(dat)){
    
    # delete the rownames of data or the PGSEA will go wrong.
    row.names(dat)<-NULL
    dat <- data.matrix(dat)
    pdata <- as.data.frame(state)
    rownames(pdata) <- colnames(dat)
    
    metadata <- data.frame(labelDescription=c("state"), row.names=c("state"))
    phenoData <- new("AnnotatedDataFrame", data=pdata, varMetadata=metadata) 
    dataExp <- ExpressionSet(assayData=dat,  phenoData=phenoData)
    dataExp
}


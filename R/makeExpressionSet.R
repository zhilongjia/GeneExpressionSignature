#' Convert data.frame into ExpressionSet Object.
#' 
#' Convert data.frame into ExpressionSet Object.
#' @param dat a data.frame with gene in row, sample in column
#' @return an ExpressionSet object
#' @author Zhilong Jia
#' @import Biobase
#' @export
#' @examples 
#' dat <- data.frame(d1=1:6, d2=6:1)
#' datExp <- makeExpressionSet(dat)
#' 

makeExpressionSet <- function(dat){
    
    # delete the rownames of data or the PGSEA will go wrong.
    row.names(dat)<-NULL
    dat <- data.matrix(dat)
    pdata <- colnames(dat)
    pdata <- as.data.frame(pdata)
    names(pdata) <- "state"
    rownames(pdata) <- pdata$state
    metadata <- data.frame(labelDescription=c("state"), row.names=c("state"))
    phenoData <- new("AnnotatedDataFrame", data=pdata, varMetadata=metadata) 
    dataExp <- ExpressionSet(assayData=dat,  phenoData=phenoData)
    dataExp
}


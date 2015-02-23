makeExpressionSet <- function(data){
    
    #delete the rownames of data or the PGSEA will go wrong.
    row.names(data)<-NULL
    data <- data.matrix(data)
    pdata <- colnames(data)
    pdata <- as.data.frame(pdata)
    names(pdata) <- "state"
    rownames(pdata) <- pdata$state
    metadata <- data.frame(labelDescription=c("state"), row.names=c("state"))
    phenoData <- new("AnnotatedDataFrame", data=pdata, varMetadata=metadata) 
    dataExp <- ExpressionSet(assayData=data,  phenoData=phenoData)
    dataExp
}


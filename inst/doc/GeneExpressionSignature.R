### R code from vignette source 'GeneExpressionSignature.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: GeneExpressionSignature.Rnw:36-37
###################################################
library(GeneExpressionSignature)


###################################################
### code chunk number 2: GeneExpressionSignature.Rnw:48-62
###################################################
# If you have network access
#GSM118720 <- getGEO('GSM118720')
# GSM118721 <- getGEO('GSM118721')
if (require(GEOquery)){
  #treatment gene-expression profiles
  GSM118720 <- getGEO(filename=system.file("extdata/GSM118720.soft",package=
  "GeneExpressionSignature"))
  #control gene-expression profiles
  GSM118721 <- getGEO(filename=system.file("extdata/GSM118721.soft",package=
  "GeneExpressionSignature"))
  #data ranking according to the different expression values 
  control <- as.matrix(as.numeric(Table(GSM118721)[,2]))
  treatment <- as.matrix(as.numeric(Table(GSM118720)[,2]))
  ranked_list <-getRLs(control,treatment)}


###################################################
### code chunk number 3: GeneExpressionSignature.Rnw:73-77
###################################################
data(exampleSet)
show(exampleSet)
exprs(exampleSet)[c(1:10),c(1:3)]
levels(as(phenoData(exampleSet),"data.frame")[,1])


###################################################
### code chunk number 4: RankMerging
###################################################
MergingSet <- RankMerging(exampleSet,"Spearman",weighted=TRUE)
show(MergingSet)


###################################################
### code chunk number 5: ScoreGSEA
###################################################
ds <- ScoreGSEA(MergingSet,250,"avg")
ds[1:5,1:5]


###################################################
### code chunk number 6: SignatureDistance
###################################################
SignatureDistance(exampleSet,SignatureLength=250,MergingDistance="Spearman",ScoringMethod="GSEA",ScoringDistance="avg",weighted=TRUE)
ds[1:5,1:5]


###################################################
### code chunk number 7: RankMerging
###################################################
MergingSet <- RankMerging(exampleSet,"Spearman",weighted=TRUE)
show(MergingSet)


###################################################
### code chunk number 8: ScoreGSEA
###################################################
ds <- ScoreGSEA(MergingSet,250,"avg")
ds[1:5,1:5]


###################################################
### code chunk number 9: ScorePGSEA
###################################################
ds <- ScorePGSEA(MergingSet,250,"avg")
ds[1:5,1:5]


###################################################
### code chunk number 10: GeneExpressionSignature.Rnw:147-152
###################################################
if (require(apcluster)){
  library(apcluster)
  clusterResult <- apcluster(1-ds)
  show(clusterResult)
}


###################################################
### code chunk number 11: GeneExpressionSignature.Rnw:163-164
###################################################
sessionInfo()



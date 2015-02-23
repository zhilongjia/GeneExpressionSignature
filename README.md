# GeneExpressionSignature

Gene Expression Signature based Similarity Metric

This package gives the implementations of the gene expression signature and its distance to each. Gene expression signature is represented as a list of genes whose expression is correlated with a biological state of interest. And its distance is defined using a nonparametric, rank-based pattern-matching strategy based on the Kolmogorov-Smirnov statistic. Gene expression signature and its distance can be used to detect similarities among the signatures of drugs, diseases, and biological states of interest.


This version (1.13.0) is based on [GeneExpressionSignature version 1.12.0; original author: Yang Cao](http://bioconductor.org/packages/release/bioc/html/GeneExpressionSignature.html)

## New functions:

**makeExpressionSet** : making ExpressionSet object from a matrix or data.frame with colnames.

*ScoreGSEA2* and *integratePRL2* are discarded due to the birth of makeExpressionSet.

## Example:

data <- makeExpressionSet(data)

ds=ScoreGSEA(MergingSet,250,"avg")


## Reference:
Cao Y (2012). GeneExpressionSignature: Gene Expression Signature based Similarity Metric. R package version 1.12.0.

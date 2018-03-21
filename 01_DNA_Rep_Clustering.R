#######################################
######### Libraries ############
#######################################
library(DESeq)
#library(arrayQualityMetrics)
library(RColorBrewer)
library(gplots)
#######################################
######### Helper Functions ############
#######################################
# Use Deseq instead
#standardize <- function(inp.vec){
#  return((inp.vec - mean(inp.vec))/sd(inp.vec))
#}

#######################################
######### Analysis ############
#######################################
dnarep.seq <- read.csv("dnarep_mat.csv", header=TRUE, row.names = 1)
dnarep.mat <- data.matrix(dnarep.seq)
#Open Question: why aren't counts integers?
dnarep.seq.rounded <- round(dnarep.mat)
#Create a "Count Data Set" object which is what the DeSeq package expects
dnarep.cds <- newCountDataSet(dnarep.seq.rounded, 
                              condition = factor(rep("Untreated", ncol(dnarep.mat))))
dnarep.cds <- estimateSizeFactors(dnarep.cds)

#Create matrix of the normalized counts
dnarep.normd <- counts(dnarep.cds, normalized=TRUE)
# Cluster the normalized data - average
dnarep.clust <- hclust(dist(dnarep.normd, method="euclidean"), 
                       method="average")
# Cluster the normalized data - complete
dnarep.clust.comp <- hclust(dist(dnarep.normd, method="euclidean"), 
                       method="complete")
# Get a matrix of variance stabilized data
dnarep.cds <- estimateDispersions( dnarep.cds, method="blind" )
dna.vsd <- getVarianceStabilizedData( dnarep.cds )
# Cluster the variance stabilized data
dnarep.vsd.clust <- hclust(dist(dna.vsd, method="euclidean"), 
                       method="average")
dnarep.vsd.clust.comp <- hclust(dist(dna.vsd, method="euclidean"),
                                method = "complete")
####
# TODO: Figure out which clustering method is most appropriate
# TODO: Figure out if we should cluster based on the normalized or variance stabilized data

# Heatmap from Deseq vignette
# This is just for trying to understand a bit more about the package
cdsFullBlind = estimateDispersions( dnarep.cds, method = "blind" )
vsdFull = varianceStabilizingTransformation( cdsFullBlind )


select = order(rowMeans(counts(dnarep.cds)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
# Cluster with the transformed data. Only takes the 3 most expressed genes
heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
# Cluster with the untransformed count data - Don't use this, just for my own edification
heatmap.2(counts(dnarep.cds)[select,], col = hmcol, trace="none", margin=c(10,6))
# Cluster with the transformed data but don't cluster samples
heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6), 
          Colv=FALSE,
          dendrogram="row",
          labCol = NULL)
# Outlier detection
#cdsBlind = estimateDispersions( dnarep.cds, method="blind" )
#vsd = varianceStabilizingTransformation( cdsBlind )
#arrayQualityMetrics(vsd)

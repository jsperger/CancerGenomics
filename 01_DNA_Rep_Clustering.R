####
# TODO: Pull in clinical data for specimen collection method
# TODO: Figure out which normalization method is most appropriate e.g. median,
#       75th quantile, one of the options from DESeq2
# TODO: Determine appropriate options for Clustering (complete vs. average), euclidean distance
# TODO: Run k-means clustering and compare to hierarchical clustering

#######################################
######### Libraries ############
#######################################
library(DESeq2)
library(ConsensusClusterPlus)
#library(arrayQualityMetrics)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library("BiocParallel")

#######################################
######### Data and Normalization ######
#######################################
###### Data input
dnarep.seq <- read.csv("dnarep_mat.csv", header=TRUE, row.names = 1)
dnarep.mat <- data.matrix(dnarep.seq)
#TCGA Barcode
# https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode
barcode <- matrix(unlist(strsplit(colnames(dnarep.mat), ".", fixed=TRUE)), ncol=7, byrow = TRUE)
colnames(barcode) <- c("Project", "TSS", "Participant", "Sample", "Portion", "Plate", "Center")
#barcode[,c(2,3,7)] <- as.numeric(barcode[,c(2,3,7)])
barcode <- data.frame(barcode)
barcode$Participant <- factor(barcode$Participant)
# Create DeSeq dataset
dna.desq <- DESeqDataSetFromMatrix(countData = dnarep.mat,
                                   colData = barcode,
                                   design = ~ 1)

#Preprocessing
# Remove genes with <10 reads, Only removes 1 gene
to.keep <- rowSums(counts(dna.desq)) >= 10
dna.desq <- dna.desq[to.keep,]
#dna.desq <- estimateSizeFactors(dna.desq)
#dna.desq <- estimateDispersions(dna.desq)
dna.desq <- DESeq(dna.desq, parallel = TRUE, BPPARAM=MulticoreParam(4))
ddr.res <- results(dna.desq)

# Basic approaches to normalization
# Median Normalization: New = Obs - median
dnarep.median.normd <- sweep(dnarep.mat,1, apply(dnarep.mat,1,median,na.rm=T))
# 75th Quantile Normalization: New = Obs - 75th quantile
dnarep.75q.normd <- sweep(dnarep.mat,1, apply(dnarep.mat,1,quantile, probs=.75,na.rm=T))

# DESeq2 Based Normalizations
# Normalizes the data according to the regularized logarithm transformation
# Note: broke R on my laptop. Will have to try running overnight or on the cluster. 
# ddr.rlog <- rlog(dna.desq)

ddr.vsd <- varianceStabilizingTransformation(dna.desq, blind = FALSE)
ddr.norm <- normTransform(dna.desq)
#######################################
######### Clustering ######
#######################################

######### Consensus Clustering
##
# Consensus Clustering based on median normalized data
# Hierarchical clustering using pearson correlation
dnarep.med.ccres <- ConsensusClusterPlus(dnarep.median.normd,maxK=6,reps=50,pItem=0.8,
                                        pFeature=1, 
                                        title="./ConsensusClustering/MedianNorm",
                                        clusterAlg="hc",
                                        distance="pearson",seed=1262118388.71279,plot="pdf")

# Consensus Clustering based on 75th Quantile normalized data
# Hierarchical clustering using pearson correlation
dnarep.med.ccres <- ConsensusClusterPlus(dnarep.75q.normd,maxK=6,reps=50,pItem=0.8,
                                         pFeature=1, 
                                         title="./ConsensusClustering/75qNorm",clusterAlg="hc",
                                         distance="pearson",seed=1262118388.71279,plot="pdf")
# Consensus Clustering based on median normalized data
# k-means clustering using euclidean distance
dnarep.med.kmeans <- ConsensusClusterPlus(dnarep.median.normd,maxK=6,reps=50,pItem=0.8,
                                         pFeature=1, 
                                         title="./ConsensusClustering/MedianNorm/kmeans",
                                         clusterAlg="km",
                                         distance="euclidean",seed=1262118388.71279,plot="pdf")
# Consensus Clustering based on 75th Quantile normalized data
# k-means clustering using euclidean distance
dnarep.75q.kmeans <- ConsensusClusterPlus(dnarep.75q.normd,maxK=6,reps=50,pItem=0.8,
                                          pFeature=1, 
                                          title="./ConsensusClustering/75qNorm/kmeans",
                                          clusterAlg="km",
                                          distance="euclidean",seed=1262118388.71279,plot="pdf")
# Consensus Clustering based on median normalized data
# Hierarchical clustering using pearson correlation
dnarep.vst.ccres <- ConsensusClusterPlus(assay(vsd),maxK=6,reps=50,pItem=0.8,
                                         pFeature=1, 
                                         title="./ConsensusClustering/VST",
                                         clusterAlg="hc",
                                         distance="euclidean",seed=1262118388.71279,plot="pdf")

# Item Consensus Plots
icl.med <- calcICL(dnarep.med.ccres,title="./ConsensusClustering/MedianNorm",plot="pdf")
icl.75q <- calcICL(dnarep.med.ccres,title="./ConsensusClustering/75qNorm",plot="pdf")


######### Hierarchical Clustering
##
# Cluster the DeSeq normalized data - average
dnarep.clust <- hclust(dist(t(dnarep.normd), method="euclidean"), 
                       method="average")
# Cluster the normalized data - complete
dnarep.clust.comp <- hclust(dist(t(dnarep.normd), method="euclidean"), 
                       method="complete")
dnarep.clust.genes <- hclust(dist(dnarep.normd, method="euclidean"), 
                                method="complete")
# TODO: Update this to DESeq2, code for DESeq not compatible
# Get a matrix of variance stabilized data
#dnarep.cds <- estimateDispersions( dnarep.cds, method="blind" )
#dna.vsd <- getVarianceStabilizedData( dnarep.cds )
# Cluster the variance stabilized data
#dnarep.vsd.clust <- hclust(dist(dna.vsd, method="euclidean"), 
#                       method="average")
#dnarep.vsd.clust.comp <- hclust(dist(dna.vsd, method="euclidean"),
#                                method = "complete")

# Heatmap from Deseq vignette
# This is just for trying to understand a bit more about the package
#cdsFullBlind = estimateDispersions( dnarep.cds, method = "blind" )
#vsdFull = varianceStabilizingTransformation( cdsFullBlind )


#select = order(rowMeans(counts(dnarep.cds)), decreasing=TRUE)[1:30]
#hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
# Cluster with the transformed data. Only takes the 3 most expressed genes
#heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
# Cluster with the untransformed count data - Don't use this, just for my own edification
#heatmap.2(counts(dnarep.cds)[select,], col = hmcol, trace="none", margin=c(10,6))
# Cluster with the transformed data but don't cluster samples
#heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6), 
#          Colv=FALSE,
#          dendrogram="row",
#          labCol = NULL)
# Cluster with all the data
#heatmap.2(exprs(vsdFull)[,], col = hmcol, trace="none")

# Outlier detection
#cdsBlind = estimateDispersions( dnarep.cds, method="blind" )
#vsd = varianceStabilizingTransformation( cdsBlind )
#arrayQualityMetrics(vsd)
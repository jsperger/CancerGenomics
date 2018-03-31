####
# TODO: Determine appropriate options for Clustering (complete vs. average), euclidean distance
# TODO: Determine if we should use a subset of genes e.g. most highly expressed

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
ddr.seq <- read.csv("ddr_mat.csv", header=TRUE, row.names = 1)
ddr.mat <- data.matrix(ddr.seq)
#TCGA Barcode
# https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode
barcode <- matrix(unlist(strsplit(colnames(ddr.mat), ".", fixed=TRUE)), ncol=7, byrow = TRUE)
colnames(barcode) <- c("Project", "TSS", "Participant", "Sample", "Portion", "Plate", "Center")
#barcode[,c(2,3,7)] <- as.numeric(barcode[,c(2,3,7)])
barcode <- data.frame(barcode)
barcode$Participant <- factor(barcode$Participant)
# Create DeSeq dataset
ddr.desq <- DESeqDataSetFromMatrix(countData = ddr.mat,
                                   colData = barcode,
                                   design = ~ TSS)

#Preprocessing
# Remove genes with <10 reads, Only removes 1 gene
to.keep <- rowSums(counts(ddr.desq)) >= 10
ddr.desq <- ddr.desq[to.keep,]
#ddr.desq <- estimateSizeFactors(ddr.desq)
#ddr.desq <- estimateDispersions(ddr.desq)
ddr.desq <- DESeq(ddr.desq, parallel = TRUE, BPPARAM=MulticoreParam(4))
ddr.res <- results(ddr.desq)

# DESeq2 Based Normalizations
# Normalizes the data according to the regularized logarithm transformation
# Note: broke R on my laptop. Will have to try running overnight or on the cluster. 
# ddr.rlog <- rlog(ddr.desq)

ddr.vsd <- varianceStabilizingTransformation(ddr.desq, blind = TRUE)
vsd=assay(ddr.vsd)  # vsd is now the normalized log2-transformed data

#######################################
######### Heat Maps / EDA ######
#######################################

# Complete sample/genelist hierarchical clustering
full.sample.cor.hmap <- pheatmap(cor(vsd), show_rownames = FALSE, show_colnames = FALSE, 
                                 main = "Sample Pearson Correlation Clustering") 
full.sample.dist.hmap <- pheatmap(dist(t(vsd)), show_rownames = FALSE, show_colnames = FALSE, 
         main = "Sample Euclidean Distance Clustering") 
# Try using spearman correlation instead of Perason
# Results look about the same
pheatmap(cor(vsd, method="spearman"), show_rownames = FALSE, show_colnames = FALSE,
         main = "Sample Spearman Correlation Clustering")  
#Cluster Method doesn't seem to be working
# These all produce idential clusterings
pheatmap(cor(vsd), show_rownames = FALSE, show_colnames = FALSE,
         cluster_method="single",
         main = "Sample Pearson Correlation - Single Linkage") 
pheatmap(cor(vsd), show_rownames = FALSE, show_colnames = FALSE, 
         clustering_method = "complete",
         main = "Sample Pearson Correlation - Complete Linkage") 
pheatmap(cor(vsd), show_rownames = FALSE, show_colnames = FALSE, 
         cluster_method="average",
         main = "Sample Pearson Correlation - Average Linkage") 

#Cluster genes using hierarchical clustering
# Using the transpose does gene correlations
pheatmap(cor(t(vsd)), show_rownames = FALSE, show_colnames = FALSE, 
         main = "Sample Pearson Correlation Transpose")  

#######################################
######### Cluster genes and samples ######
#######################################

# Function for selecting genes and creating a heatmap which clusters rows and columns
# Methods: most expressed, most varying, highest coefficient of variation, both expressed and varying
# n.genes sets the number of genes you want to select 
#  except for combined where the number of genes will be <= 2*n.genes
# Changes some of the defaults of pheatmap
phs <- function(inp.desq = ddr.vsd, smethod = "mean", n.genes = 30, ...){
  stopifnot(any(smethod == c("mean", "var", "cv", "combined", "intersect")))
  ### Different methods for selecting genes
  #Selects the most expressed genes
  if(smethod == "mean"){
    select <- order(rowMeans(assay(inp.desq)),
                         decreasing=TRUE)[1:n.genes]
  }
  # Selects the genes with the highest variance
  # Not sure the trade off between sample size and wanting genes which actually vary
  if(smethod == "var"){
    select <- order(rowVars(assay(inp.desq)),
                    decreasing=TRUE)[1:n.genes]
  }
  # Selects the genes with the highest coefficient of variation
  # cv = variance / mean
  if(smethod == "cv"){
    select <- order(rowVars(assay(inp.desq))/rowMeans(assay(inp.desq)),
                    decreasing=TRUE)[1:n.genes]
  }
  # Selects both high variance and highly expressed genes
  if(smethod == "combined"){
    mean.select <- order(rowMeans(assay(inp.desq)),
                    decreasing=TRUE)[1:n.genes]
    var.select <- order(rowVars(assay(inp.desq)),
                    decreasing=TRUE)[1:n.genes]
    select <- unique(c(mean.select, var.select))
  }
  # Selects both high variance and highly expressed genes
  if(smethod == "intersect"){
    mean.select <- order(rowMeans(assay(inp.desq)),
                         decreasing=TRUE)[1:n.genes]
    var.select <- order(rowVars(assay(inp.desq)),
                        decreasing=TRUE)[1:n.genes]
    select <- mean.select[mean.select %in% var.select]
  }
  
  pheatmap(assay(inp.desq)[select,], ...)
}

#Heatmap and hierarchical clustering of 40 most expressed genes
phs(smethod="mean", n.genes = 40, main="Clustering 40 Most Expressed Genes", 
    cluster_rows=TRUE, show_rownames=TRUE, 
    cluster_cols=TRUE, show_colnames=FALSE)
#Heatmap and hierarchical clustering of 40 most varying genes
phs(smethod="var", n.genes = 40, main="Clustering 40 Most Varying Genes", 
    cluster_rows=TRUE, show_rownames=TRUE, 
    cluster_cols=TRUE, show_colnames=FALSE)

#Heatmap and hierarchical clustering of 40 most expressed and most varying genes
phs(smethod="combined", n.genes = 40, main="Clustering 40 Most Expressed and Varying Genes", 
    cluster_rows=TRUE, show_rownames=TRUE, 
    cluster_cols=TRUE, show_colnames=FALSE)

#Heatmap and hierarchical clustering of 50 genes with highest cv
# Cluster for genes is based on correlation
# Cluster for samples is based on euclidean distance
phs(smethod="cv", n.genes = 50, main="50 Genes with Highest CV - Genes Cor, Samples Distance", 
    cluster_rows=TRUE, show_rownames=TRUE, 
    cluster_cols=TRUE, show_colnames=FALSE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "euclidean")
#Heatmap and hierarchical clustering of 50 most expressed and most varying genes
# Cluster for genes and samples is based on euclidean distance
phs(smethod="cv", n.genes = 50, main="50 Genes with Highest CV - Genes and Samples Distance", 
    cluster_rows=TRUE, show_rownames=TRUE, 
    cluster_cols=TRUE, show_colnames=FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean")

#######################################
######### Cluster Assignments ######
#######################################
#Compute the euclidean distance of the vst-transformed samples
vsd.dist <- dist(t(vsd), method = "euclidean")
#Hierarchical clustering of the samples using the single linkage method
hc.single <- hclust(vsd.dist, method="single")
#Hierarchical clustering of the samples using the average linkage method
hc.avg <- hclust(vsd.dist, method="average")
#Hierarchical clustering of the samples using the complete linkage method
hc.comp <- hclust(vsd.dist, method="complete")
#Hierarchical clustering of the samples using the centroid  method
hc.cent <- hclust(vsd.dist, method="centroid")
# TODO: Also want to look at what happens if you try and cluster with a smaller number of genes
# How much changes
#hc.cvs.comp <- hclust(dist(t(vsd[cv.select,])), method="complete")
#hc.cvs.comp.asgn <- sapply(2:10, cutree, tree=hc.cvs.comp, h=NULL)
# Let's see how the cluster sizes vary across methods
# Only complete assignment creates clusters with more than a couple members in each
# That's not to say it's the correct method
hc.single.asgn <- sapply(2:10, cutree, tree=hc.single, h=NULL)
hc.avg.asgn <- sapply(2:10, cutree, tree=hc.avg, h=NULL)
hc.comp.asgn <- sapply(2:10, cutree, tree=hc.comp, h=NULL)
hc.cent.asgn <- sapply(2:10, cutree, tree=hc.cent, h=NULL)
#apply(hc.cent.asgn, 2, table)
#######################################
######### Consensus Clustering ######
#########################################

# Consensus Clustering based on VST normalized data
# Hierarchical clustering using pearson correlation
ddr.vst.ccres <- ConsensusClusterPlus(vsd,maxK=10,reps=50,pItem=0.8,
                                         pFeature=1, 
                                         title="./ConsensusClustering/VST",
                                         clusterAlg="hc",
                                      innerLinkage = "complete",
                                      finalLinkage = "complete",
                                         distance="euclidean",seed=1262118388.71279,plot="pdf")

# Item Consensus Plots
icl.vst <- calcICL(ddr.vst.ccres,title="./ConsensusClustering/VST",plot="pdf")

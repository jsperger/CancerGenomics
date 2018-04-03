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
library(flexclust)
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

cvs <- sqrt(rowVars(assay(ddr.desq)))/rowMeans(assay(ddr.desq))
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
######### Cluster Heatmaps ######
#######################################

# Function for selecting genes and creating a heatmap which clusters rows and columns
# Methods: most expressed, most varying, highest coefficient of variation, both expressed and varying
# n.genes sets the number of genes you want to select 
#  except for combined where the number of genes will be <= 2*n.genes
# Changes some of the defaults of pheatmap
phs <- function(sort.desq = ddr.desq, inp.vst = ddr.vsd, smethod = "mean", n.genes = 30, ord.dec = TRUE, ...){
  stopifnot(any(smethod == c("mean", "var", "cv", "combined", "intersect")))
  ### Different methods for selecting genes
  #Selects the most expressed genes
  if(smethod == "mean"){
    select <- order(rowMeans(assay(sort.desq)),
                         decreasing=ord.dec)[1:n.genes]
  }
  # Selects the genes with the highest variance
  # Not sure the trade off between sample size and wanting genes which actually vary
  if(smethod == "var"){
    select <- order(rowVars(assay(sort.desq)),
                    decreasing=ord.dec)[1:n.genes]
  }
  # Selects the genes with the highest coefficient of variation
  # cv = variance / mean
  if(smethod == "cv"){
    select <- order(sqrt(rowVars(assay(sort.desq)))/rowMeans(assay(sort.desq)),
                    decreasing=ord.dec)[1:n.genes]
  }
  # Selects both high variance and highly expressed genes
  if(smethod == "combined"){
    mean.select <- order(rowMeans(assay(sort.desq)),
                    decreasing=ord.dec)[1:n.genes]
    var.select <- order(rowVars(assay(sort.desq)),
                    decreasing=ord.dec)[1:n.genes]
    select <- unique(c(mean.select, var.select))
  }
  # Selects both high variance and highly expressed genes
  if(smethod == "intersect"){
    mean.select <- order(rowMeans(assay(sort.desq)),
                         decreasing=ord.dec)[1:n.genes]
    var.select <- order(rowVars(assay(sort.desq)),
                        decreasing=ord.dec)[1:n.genes]
    select <- mean.select[mean.select %in% var.select]
  }
  
  pheatmap(assay(inp.vst)[select,], ...)
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
phs(smethod="cv", n.genes = 50, main="50 Genes with Highest Coef of Variation sigam/mu", 
    cluster_rows=TRUE, show_rownames=TRUE, 
    cluster_cols=TRUE, show_colnames=FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean")
phs(smethod="cv", n.genes = 50, main="50 Genes with Lowest Coef of Variation sigam/mu", 
    cluster_rows=TRUE, show_rownames=TRUE, 
    cluster_cols=TRUE, show_colnames=FALSE,
    ord.dec = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean")
phs(smethod="mean", n.genes = 200, main="200 Genes with Highest Expression", 
    cluster_rows=TRUE, show_rownames=FALSE, 
    cluster_cols=TRUE, show_colnames=FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean")
phs(smethod="cv", n.genes=50, main = "50 Genes with ")

test.cv <- order(sqrt(rowVars(assay(ddr.desq)))/rowMeans(assay(ddr.desq)),
                 decreasing=TRUE)[1:50]
alt.cv <- order(sqrt(rowVars(assay(ddr.desq)))/rowMeans(assay(ddr.desq)),
                   decreasing=FALSE)[1:50]
pheatmap(vsd[test.cv,],main="50 Genes with Highest CV from original assay - decreasing", 
         cluster_rows=TRUE, show_rownames=TRUE, 
         cluster_cols=TRUE, show_colnames=FALSE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

pheatmap(vsd[alt.cv,],main="50 Genes with Highest CV from original assay - increasing", 
         cluster_rows=TRUE, show_rownames=TRUE, 
         cluster_cols=TRUE, show_colnames=FALSE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

test.mean1 <- order(rowMeans(assay(ddr.desq)),decreasing=TRUE)[1:40]
test.mean2 <- order(rowMeans(assay(ddr.vsd)),decreasing=TRUE)[1:40]
pheatmap(vsd[test.mean1,],main="50 Genes with Expression from original assay", 
         clustering_method = "ward.D",
         cluster_rows=TRUE, show_rownames=TRUE, 
         cluster_cols=TRUE, show_colnames=FALSE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")
pheatmap(vsd[test.mean2,],main="40 Genes with Highest Mean Expression", 
         cluster_rows=TRUE, show_rownames=TRUE, 
         cluster_cols=TRUE, show_colnames=FALSE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")
test.mean3 <- order(rowMeans(assay(ddr.desq)),decreasing=TRUE)[1:100]
pheatmap(vsd[test.mean3,],main="100 Genes with Expression from original assay", 
         cluster_rows=TRUE, show_rownames=TRUE, 
         cluster_cols=TRUE, show_colnames=FALSE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

test.cv.forint <- order(sqrt(rowVars(assay(ddr.desq)))/rowMeans(assay(ddr.desq)),
                        decreasing=FALSE)[1:100]
mean.cv.int <- test.mean3[test.mean3 %in% test.cv.forint]
pheatmap(vsd[mean.cv.int,],main="Genes with Low CV and High Expression", 
         cluster_rows=TRUE, show_rownames=TRUE, 
         cluster_cols=TRUE, show_colnames=FALSE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")
ts <- ddr.res$padj < .1 & !is.na(ddr.res$padj)
difexp <- ddr.vsd[ts,]
pheatmap(vsd[ts,],main="Sig Differentially expressed", 
         cluster_rows=TRUE, show_rownames=TRUE, 
         cluster_cols=TRUE, show_colnames=FALSE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")
#######################################
######### Cluster Methods ######
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

hc.mean.comp <- hclust(dist(t(vsd[test.mean1,])), method="complete")
hc.mean.avg <- hclust(dist(t(vsd[test.mean1,])), method="average")
hc.mean.cent <- hclust(dist(t(vsd[test.mean1,])), method="centroid")
hc.mean.ward <- hclust(dist(t(vsd[test.mean1,])), method="ward.D")
km <- kmeans(t(vsd[test.mean1,]), centers=3, iter.max=30, nstart=50)
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
hc.mean.asgn <- sapply(2:10, cutree, tree=hc.mean.comp, h=NULL)
hc.mean.cent.asgn <- sapply(2:10, cutree, tree=hc.mean.cent, h=NULL)
hc.mean.avg.asgn <- sapply(2:10, cutree, tree=hc.mean.avg, h=NULL)
hc.mean.ward.asgn <- sapply(2:10, cutree, tree=hc.mean.ward, h=NULL)


apply(hc.comp.asgn, 2, table)
apply(hc.mean.cent.asgn, 2, table)
apply(hc.mean.avg.asgn, 2, table)
apply(hc.mean.ward.asgn, 2, table)

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
ddr.vst.ccres <- ConsensusClusterPlus(vsd[test.mean1,],maxK=10,reps=50,pItem=0.8,
                                      pFeature=1, 
                                      title="./ConsensusClustering/VST/Mean",
                                      clusterAlg="hc",
                                      innerLinkage = "complete",
                                      finalLinkage = "complete",
                                      distance="euclidean",seed=1262118388.71279,plot="pdf")
ddr.vst.ccres <- ConsensusClusterPlus(vsd[test.mean1,],maxK=10,reps=50,pItem=0.8,
                                      pFeature=1, 
                                      title="./ConsensusClustering/VST/Mean/Ward",
                                      clusterAlg="hc",
                                      innerLinkage = "ward.D",
                                      finalLinkage = "ward.D",
                                      distance="euclidean",seed=1262118388.71279,plot="pdf")
ConsensusClusterPlus(vsd,maxK=10,reps=50,pItem=0.8,
                     pFeature=1, 
                     title="./ConsensusClustering/VST/Pearson/Complete",
                     clusterAlg="hc",
                     innerLinkage = "complete",
                     finalLinkage = "complete",
                     distance="pearson",seed=1262118388.71279,plot="pdf")
ConsensusClusterPlus(vsd,maxK=10,reps=50,pItem=0.8,
                     pFeature=1, 
                     title="./ConsensusClustering/VST/Pearson/Complete",
                     clusterAlg="hc",
                     innerLinkage = "complete",
                     finalLinkage = "complete",
                     distance="pearson",seed=1262118388.71279,plot="pdf")

######### k-Means Consensus Clustering ######
kmpp <- kcca(vsd[test.mean1,], k=3, family=kccaFamily("kmeans"),
           control=list(initcent="kmeanspp"), save.data = TRUE)
kmeans.cc <- ConsensusClusterPlus(vsd[test.mean1,],maxK=10,reps=50,pItem=0.8,
                                      pFeature=1, 
                                      title="./ConsensusClustering/VST/kMeans",
                                      clusterAlg="km",
                                      distance="euclidean",seed=1262118388.71279,plot="pdf")

difexp.cc <- ConsensusClusterPlus(vsd[ts,],maxK=10,reps=50,pItem=0.8,
                                  pFeature=1, 
                                  title="./ConsensusClustering/VST/Difexp",
                                  clusterAlg="km",
                                  distance="euclidean",seed=1262118388.71279,plot="pdf")
difexp.hc.cc <- ConsensusClusterPlus(vsd[ts,],maxK=10,reps=50,pItem=0.8,
                                  pFeature=1, 
                                  title="./ConsensusClustering/VST/Difexp/HC",
                                  clusterAlg="hc",
                                  innerLinkage = "complete",
                                  finalLinkage = "complete",
                                  distance="euclidean",seed=1262118388.71279,plot="pdf")
difexp.pam.cc <- ConsensusClusterPlus(vsd[ts,],maxK=10,reps=50,pItem=0.8,
                                     pFeature=1, 
                                     title="./ConsensusClustering/VST/Difexp/pam",
                                     clusterAlg="pam",
                                     distance="euclidean",seed=1262118388.71279,plot="pdf")
# Item Consensus Plots
icl.vst <- calcICL(ddr.vst.ccres,title="./ConsensusClustering/VST",plot="pdf")

#######################################
######### Cluster Assignments ######
#######################################
cluster.assignments <- data.frame(Patient = paste0(barcode$Project,"-",barcode$TSS,"-",barcode$Participant),
                      hcC2 = hc.mean.asgn[,1],
                      hcC3 = hc.mean.asgn[,2],
                      hcC4 = hc.mean.asgn[,3])

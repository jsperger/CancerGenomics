####
# TODO: Pull in clinical data for specimen collection method
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
######### Clustering ######
#######################################

#Heatmap and hierarchical clustering of 40 most differentially expressed genes
mean.select <- order(rowMeans(assay(ddr.vsd)),
                decreasing=TRUE)[1:40]
pheatmap(assay(ddr.vsd)[mean.select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, show_colnames = FALSE)

var.select <- order(rowVars(assay(ddr.vsd)),
                     decreasing=TRUE)[1:40]
pheatmap(assay(ddr.vsd)[var.select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, show_colnames = FALSE)

combined.select <- unique(c(mean.select, var.select))
inter.select <- mean.select[mean.select %in%  var.select]
pheatmap(assay(ddr.vsd)[inter.select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, show_colnames = FALSE, 
         main="High Expression High Variance - Euclidean")
pheatmap(assay(ddr.vsd)[inter.select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, show_colnames = FALSE, 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main="High Expression High Variance - Correlation")
full.sample.cor.hmap <- pheatmap(cor(vsd), show_rownames = FALSE, show_colnames = FALSE, 
         main = "Sample Pearson Correlation")  
pheatmap(cor(vsd, method="spearman"), show_rownames = FALSE, show_colnames = FALSE,
         main = "Sample Spearman Correlation")  
coef.var <- rowVars(assay(ddr.vsd))/rowMeans(assay(ddr.vsd))
cv.select <- order(coef.var,
                    decreasing=TRUE)[1:300]
pheatmap(cor(assay(ddr.vsd)[cv.select,]), 
         show_rownames=FALSE, show_colnames = FALSE, 
         main="Coefficient of Variation Selected Heatmap")
#######################################
######### Consensus Clustering ######
#########################################

# Consensus Clustering based on VST normalized data
# Hierarchical clustering using pearson correlation
ddr.vst.ccres <- ConsensusClusterPlus(vsd,maxK=10,reps=50,pItem=0.8,
                                         pFeature=1, 
                                         title="./ConsensusClustering/VST",
                                         clusterAlg="hc",
                                         distance="euclidean",seed=1262118388.71279,plot="pdf")
ddr.vst.ccres2 <- ConsensusClusterPlus(vsd,maxK=10,reps=50,pItem=0.8,
                                      pFeature=1, 
                                      title="./ConsensusClustering/VST/Pearson",
                                      clusterAlg="hc",
                                      distance="pearson",seed=1262118388.71279,plot="pdf")
# Item Consensus Plots
icl.vst <- calcICL(ddr.vst.ccres,title="./ConsensusClustering/VST",plot="pdf")
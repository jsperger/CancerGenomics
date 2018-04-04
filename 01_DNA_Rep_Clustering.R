#######################################
######### Libraries ############
#######################################
library("SummarizedExperiment")
library(DESeq2)
library(NMF)
library("BiocParallel")
library(pheatmap)
#library(CancerSubtypes)
#library(ConsensusClusterPlus)
#library(RColorBrewer)
#library(gplots)
#######################################
######### Data and Normalization ######
#######################################
###### Data input
ddr.seq <- read.csv("ddr_mat.csv", header=TRUE, row.names = 1)
ddr.mat <- data.matrix(ddr.seq)
#Preprocessing
# Remove genes with <10 reads, Only removes 1 gene
reads.geq10 <- (apply(ddr.mat,1,sum) >= 10)
ddr.desq <- ddr.mat[reads.geq10,]
#TCGA Barcode
# https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode
barcode <- matrix(unlist(strsplit(colnames(ddr.mat), ".", fixed=TRUE)), ncol=7, byrow = TRUE)
colnames(barcode) <- c("Project", "TSS", "Participant", "Sample", "Portion", "Plate", "Center")
barcode <- data.frame(barcode)
barcode$Participant <- factor(barcode$Participant)
#This matches the variable of the same name in the clinical data
barcode$bcr_patient_barcode = paste0(barcode$Project,"-",barcode$TSS,"-",barcode$Participant)

# DESeq2-Based Normalization: Variance Stabilizing Transformation
# ddr.vsd is now the regularized log2-transformed data
ddr.vsd <- varianceStabilizingTransformation(ddr.mat, blind = TRUE)

#######################################
######### NMF Clustering ######
#######################################

#######################################
# The following code takes a long time to run (300 runs took 1.5 hours to run on 4 cores on my laptop)
#estim.r <- nmf(ddr.vsd, 2:6, nrun=30, .opt='vp4')
#save(estim.r, file='nmf_k_comparison.RData')
#plot(estim.r)
#consensusmap(estim.r,  labCol=NA, labRow=NA)
#vsd.nmf <- nmf(x=ddr.vsd, rank=3, nrun=300, .opt='vp4')
#save(vsd.nmf, file='./Results/nmf_k3_300.RData')
load("./Results/nmf_k3_300.RData")

# Get details about the meta-genes
meta.genes <- extractFeatures(vsd.nmf)
c1.genes <- rownames(ddr.vsd)[meta.genes[[1]]]
c2.genes <- rownames(ddr.vsd)[meta.genes[[2]]]
c3.genes <- rownames(ddr.vsd)[meta.genes[[3]]]
write.table(paste0(c1.genes, collapse=", "), file="./Results/metagene_1_list.txt", 
            row.names=FALSE, col.names=FALSE, eol= "", quote = FALSE)
write.table(paste0(c2.genes, collapse=", "), file="./Results/metagene_2_list.txt", 
            row.names=FALSE, col.names=FALSE, eol= "",  quote = FALSE)
write.table(paste0(c3.genes, collapse=", "), file="./Results/metagene_3_list.txt", 
            row.names=FALSE, col.names=FALSE, eol= "",  quote = FALSE)

#######################################
######### Cluster Assignments ######
#######################################

#Assign Clusters from k=3 NMF Clustering
barcode$NMFC3 <- predict(vsd.nmf)

write.csv(barcode, file = "./Results/cluster_results.csv", row.names = FALSE)

#######################################
######### DESeq2 ######
#######################################

# Create DeSeq dataset
ddr.desq <- DESeqDataSetFromMatrix(countData = ddr.mat,
                                   colData = barcode,
                                   design = ~ NMFC3)


#Calculate Results
ddr.desq <- DESeq(ddr.desq, parallel = TRUE, BPPARAM=MulticoreParam(4))
#save(ddr.desq, file='./Results/ddr.desq.Rdata')

ddr.res <- results(ddr.desq)
summary(ddr.res, alpha=.05)

#######################################
######### Plots ######
#######################################
pdf("basis_map.pdf")
basismap(vsd.nmf, subsetRow = TRUE)
dev.off()
pdf("coef_map.pdf")
coefmap(vsd.nmf)
dev.off()
pdf("consensus_map.pdf")
consensusmap(vsd.nmf)
dev.off()
#######################################
######### Old ######
#######################################
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

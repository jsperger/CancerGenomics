#######################################
####### Existing Subtype Replication ##
#######################################

# this is an extra package that was used in the Nature paper
# install.packages('NMF')
library(NMF)

# load the Consensus Clustering function
library(ConsensusClusterPlus)

#######################################
####### Data Setup ####################
#######################################

# read the data in -- John's start --
#
# first we load the 1500 genes used for the clustering 
tcga <- data.matrix(read.table("TCGA_489_UE.top1500.txt", header=TRUE))
#
# now we'll read the data into memory (should be done in 00_Libs)
#
##
## full.rna.seq.mat <- assay(rna.seqdat)
## rownames(full.rna.seq.mat) <- rowData(rna.seqdat)$external_gene_name
## write.csv(full.rna.seq.mat, file = "full_rna_seq_mat.csv",
## row.names = TRUE)
#
# here we create an index of gene names (done in 00_Libs)
#
## gene.names <- rowData(rna.seqdat)$external_gene_name
#
# let's find the the subset of the 1500 we still have in our data
to.keep <- gene.names %in% rownames(tcga)
subtype.mat <- full.rna.seq.mat[to.keep,]
dim(subtype.mat)
# should be 1357 remaining for 379 samples

#######################################
####### Data Preprocessing ############
#######################################

# now we'll do the barcoding -- Borrowed from John --
# in the line below, I had to change "." to "-" as the strsplit delimiter
barcode <- matrix(unlist(strsplit(colnames(subtype.mat), "-", fixed=TRUE)), ncol=7, byrow = TRUE)
colnames(barcode) <- c("Project", "TSS", "Participant", "Sample", "Portion", "Plate", "Center")
#barcode[,c(2,3,7)] <- as.numeric(barcode[,c(2,3,7)])
barcode <- data.frame(barcode)
barcode$Participant <- factor(barcode$Participant)
#This matches the variable of the same name in the clinical data
barcode$bcr_patient_barcode = paste0(barcode$Project,"-",barcode$TSS,"-",barcode$Participant)

# DESeq2-Based Normalization: Variance Stabilizing Transformation
# subtype.vsd is now the regularized log2-transformed data
subtype.vsd <- varianceStabilizingTransformation(subtype.mat, blind = TRUE)

# Danger Zone

#######################################
######### NMF Clustering ######
#######################################

#######################################
# You prolly do not want to run the next couple of lines
# t <- Sys.time()
# subtype_estim.r <- nmf(subtype.vsd, 2:6, nrun=30, .opt='vp4')
# save(subtype_estim.r, file='nmf_k_subtype.RData')
# plot(subtype_estim.r)
# consensusmap(subtype_estim.r,  labCol=NA, labRow=NA)
# t2 <- Sys.time()
# t2 - t
###
t <- Sys.time()
subtype_vsd.nmf <- nmf(x=subtype.vsd, rank=4, nrun=300, .opt='vp4')
save(subtype_vsd.nmf, file='CancerGenomics/Results/nmf_k4_subtype_300.RData')
plot(subtype_vsd.nmf)
consensusmap(subtype_vsd.nmf,  labCol=NA, labRow=NA)
#load("./Results/nmf_k4_subtype_300.RData")
t2 <- Sys.time()
t2 - t

# Get details about the meta-genes
subtype_meta.genes <- extractFeatures(subtype_vsd.nmf)
sub.c1.genes <- rownames(subtype.vsd)[subtype_meta.genes[[1]]]
sub.c2.genes <- rownames(subtype.vsd)[subtype_meta.genes[[2]]]
sub.c3.genes <- rownames(subtype.vsd)[subtype_meta.genes[[3]]]
sub.c4.genes <- rownames(subtype.vsd)[subtype_meta.genes[[4]]]
write.table(paste0(sub.c1.genes, collapse=", "), file="CancerGenomics/Results/subtype_metagene_1_list.txt", 
            row.names=FALSE, col.names=FALSE, eol= "", quote = FALSE)
write.table(paste0(sub.c2.genes, collapse=", "), file="CancerGenomics/Results/subtype_metagene_2_list.txt", 
            row.names=FALSE, col.names=FALSE, eol= "",  quote = FALSE)
write.table(paste0(sub.c3.genes, collapse=", "), file="CancerGenomics/Results/subtype_metagene_3_list.txt", 
            row.names=FALSE, col.names=FALSE, eol= "",  quote = FALSE)
write.table(paste0(sub.c4.genes, collapse=", "), file="CancerGenomics/Results/subtype_metagene_4_list.txt", 
            row.names=FALSE, col.names=FALSE, eol= "",  quote = FALSE)

#######################################
######### Cluster Assignments ######
#######################################

#Assign Clusters from k=4 NMF Clustering
barcode$NMFC3 <- predict(subtype_vsd.nmf)

write.csv(barcode, file = "CancerGenomics/Results/subtype_cluster_results.csv", row.names = FALSE)

#######################################
######### DESeq2 ######
#######################################

# Create DeSeq dataset
subtype.desq <- DESeqDataSetFromMatrix(countData = subtype.mat,
                                   colData = barcode,
                                   design = ~ NMFC3)

# parallelization stuff
# source("https://bioconductor.org/biocLite.R")
# biocLite('MulticoreParam')
#library('MulticoreParam')
#Calculate Results
subtype.desq <- DESeq(subtype.desq)#, parallel = TRUE, BPPARAM=MulticoreParam(4))
#save(subtype.desq, file='CancerGenomics/Results/subtype.desq.Rdata')

subtype.res <- results(subtype.desq)
summary(subtype.res, alpha=.05)

#######################################
######### Plots ######
#######################################
pdf("sub_basis_map.pdf")
basismap(subtype_vsd.nmf, subsetRow = TRUE)
dev.off()
pdf("sub_coef_map.pdf")
coefmap(subtype_vsd.nmf)
dev.off()
pdf("sub_consensus_map.pdf")
consensusmap(subtype_vsd.nmf)
dev.off()
#######################################

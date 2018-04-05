#######################################
####### Existing Subtype Replication ##
#######################################

# this is an extra package that was used in the Nature paper - not necessary
# install.packages('NMF')
# library(NMF)

# load the Consensus Clustering function
library(ConsensusClusterPlus)

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
#
# now we'll do the barcoding -- Borrowed from John --
barcode <- matrix(unlist(strsplit(colnames(subtype.mat), ".", fixed=TRUE)), ncol=7, byrow = TRUE)
colnames(barcode) <- c("Project", "TSS", "Participant", "Sample", "Portion", "Plate", "Center")
#barcode[,c(2,3,7)] <- as.numeric(barcode[,c(2,3,7)])
barcode <- data.frame(barcode)
barcode$Participant <- factor(barcode$Participant)
# Create DeSeq dataset
subtype.desq <- DESeqDataSetFromMatrix(countData = ddr.mat,
                                   colData = barcode,
                                   design = ~ TSS)
#
# And some data cleaning -- Borrowed from John --
# Remove genes with <10 reads, Only removes 1 gene
to.keep <- rowSums(counts(subtype.desq)) >= 10
subtype.desq <- subtype.desq[to.keep,]
subtype.desq <- DESeq(ddr.subtype)#, parallel = TRUE, BPPARAM=MulticoreParam(4)) # Parallelize if you wish
subtype.res <- results(subtype.desq)

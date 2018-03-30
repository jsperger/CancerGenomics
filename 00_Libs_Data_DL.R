#######################################
######### Libraries ############
#######################################

### Install and Load Required Packages

source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
biocLite("RTCGAToolbox")
biocLite("SummarizedExperiment")
biocLite("BiocParallel")
biocLite("DESeq2")
biocLite("arrayQualityMetrics")
biocLite("ConsensusClusterPlus")
library("TCGAbiolinks")
library("RTCGAToolbox")
library("SummarizedExperiment")

#######################################
######### Data Download ############
#######################################
### Download RNA Seq Data
### Function automatically will check if the files have been downloaded
#Download Size is ~100MB
rna.query <- GDCquery(project = "TCGA-OV",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts")
GDCdownload(rna.query, method = "api", files.per.chunk = 10)

### Download Clinical Data
### Function automatically will check if the files have been downloaded
#Download Size is ~32MB
clinical.query <- GDCquery(project = "TCGA-OV", 
                           data.category = "Clinical")
GDCdownload(clinical.query)

#######################################
######### Data Preparation ############
#######################################
### Prepare RNA Data for Analysis
# Prepares a 'Summarized Experiment' object
# Writes the 'Summarized Experiment' object to an R data file "rnaSeq.rda"
rna.seqdat <- GDCprepare(query = rna.query, 
                      save = TRUE, 
                      save.filename = "rnaSeq.rda")

# Create a matrix of the RNA-seq data with genes as rows and samples as columns
full.rna.seq.mat <- assay(rna.seqdat)
rownames(full.rna.seq.mat) <- rowData(rna.seqdat)$external_gene_name

# Write the RNA-Seq data to a matrix
write.csv(full.rna.seq.mat, file = "full_rna_seq_mat.csv",
          row.names = TRUE)

### Subset the RNA-Seq data to include only genes involved in DNA repair 
# Pull the list of genes in the full data set
gene.names <- rowData(rna.seqdat)$external_gene_name

# Read in the list of genes associated with DNA Damage Repair
ddr.gene.list <- read.csv("dna_repair_genelist.csv", stringsAsFactors = FALSE)
# Take the genes which appear in the RNA-Seq data
to.keep <- gene.names %in% ddr.gene.list$Gene
ddr.mat <- full.rna.seq.mat[to.keep,]

# Write the RNA-Seq data to a matrix
write.csv(ddr.mat, file = "ddr_mat.csv",
          row.names = TRUE)
# Gets a list of which DDR Genes appear in the data
genes.in.data <- ddr.gene.list[ddr.gene.list$Gene %in% rownames(full.rna.seq.mat),]
write.csv(genes.in.data, file = "ddr.genes.in.data.csv",
          row.names = TRUE)
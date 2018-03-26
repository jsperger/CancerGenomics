#######################################
######### Libraries ############
#######################################

### Install and Load Required Packages

source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
biocLite("RTCGAToolbox")
biocLite("SummarizedExperiment")
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

### Subset the RNA-Seq data to include only genes involved in DNA repair 
# Pull the list of genes in the full data set
full.gene.list <- row.names(full.rna.seq.mat)
# The gene name format is Name|##### not sure what the numbers mean
# This splits into separate characters for Name and Number, then just takes the names
# Ugly code but I forget regular expressions
gene.list.split <- matrix(unlist(strsplit(full.gene.list, "|", fixed=TRUE)), ncol=2, byrow=TRUE)
full.gene.names <- gene.list.split[,1]
# Replace the rownames with just the gene names
rownames(full.rna.seq.mat) <- full.gene.names

# Write the RNA-Seq data to a matrix
write.csv(full.rna.seq.mat, file = "full_rna_seq_mat.csv",
          row.names = TRUE)

# Read in the list of genes associated with DNA Damage Repair
dnarep.gene.list <- read.csv("dna_repair_genelist.csv", stringsAsFactors = FALSE)
# Take the genes which appear in the RNA-Seq data
genes.in.data <- dnarep.gene.list[dnarep.gene.list$Gene %in% rownames(full.rna.seq.mat),]
dnarep.mat <- full.rna.seq.mat[genes.in.data,]

# Write the RNA-Seq data to a matrix
write.csv(dnarep.mat, file = "dnarep_mat.csv",
          row.names = TRUE)

#######################################
######### Deprecated ############
#######################################
### Download DNA Methylation Data
# Function automatically will check if the files have been downloaded
#Download Size is ~210MB

#methyl.query <- GDCquery(project = c("TCGA-OV"),
#                  data.category = "DNA methylation",
#                  platform = "Illumina Human Methylation 450",
#                  legacy = TRUE)
#GDCdownload(methyl.query, method = "api", files.per.chunk = 10)



#Loads clinical data into R. 
#TODO: Move to another script
#clinical <- GDCprepare_clinic(query, clinical.info = "patient")
#######################################
######### Libraries ############
#######################################

### Install and Load Required Packages

source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
biocLite("RTCGAToolbox")
biocLite("SummarizedExperiment")
biocLite("DESeq")
library("TCGAbiolinks")
library("RTCGAToolbox")
library("SummarizedExperiment")

#######################################
######### Data Download ############
#######################################
### Download RNA Seq Data
### Function automatically will check if the files have been downloaded
#Download Size is ~130MB
rna.query <- GDCquery(project = "TCGA-OV",
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  #platform = "Illumina HiSeq", 
                  file.type  = "results",
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE)
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
# Write the RNA-Seq data to a matrix
write.csv(full.rna.seq.mat, file = "full_rna_seq_mat.csv",
           row.names = TRUE,
           col.names = TRUE)

### Subset the RNA-Seq data to include only genes involved in DNA repair 
dnarep.gene.list <- read.csv("dna_repair_genelist.csv")
dnarep.mat <- full.rna.seq.mat[dnarep.gene.list$Gene,]
# Write the RNA-Seq data to a matrix
write.csv(dnarep.mat, file = "dnarep_mat.csv",
          row.names = TRUE,
          col.names = TRUE)

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
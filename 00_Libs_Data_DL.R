#######################################
######### Data Preparation ############
#######################################

### Install and Load Required Packages

source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
biocLite("RTCGAToolbox")
library("TCGAbiolinks")
library("RTCGAToolbox")

### Download RNA Seq Data
### Function automatically will check if the files have been downloaded
#Download Size is ~130MB
rna.query <- GDCquery(project = "TCGA-OV",
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq", 
                  file.type  = "normalized_results",
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE)
GDCdownload(rna.query, method = "api", files.per.chunk = 10)

### Download Clinical Data
### Function automatically will check if the files have been downloaded
#Download Size is ~32MB
clinical.query <- GDCquery(project = "TCGA-OV", 
                  data.category = "Clinical")
GDCdownload(clinical.query)

#Loads clinical data into R. 
#TODO: Move to another script
#clinical <- GDCprepare_clinic(query, clinical.info = "patient")
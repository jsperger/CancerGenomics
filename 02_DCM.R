#######################################
######### Libraries and Data ##########
#######################################
dcm.file.list <- list.files("./DCM/R", full.names = TRUE)
source_files(dcm.file.list)

load("./Results/nmf_k3_300.RData")
barcode <- read.csv("./Results/cluster_results.csv")
#######################################
######### Run DCM ############
#######################################
#Separate Clusters
clust1 <- ddr.vsd[,barcode$NMFC3 == 1]
clust2 <- ddr.vsd[,barcode$NMFC3 == 2]
clust3 <- ddr.vsd[,barcode$NMFC3 == 3]

# Run DCM Pairwise on Clusters
# Alpha level set using Bonferroni Correction
# TODO: Figure out if we should run with resid.full=TRUE or not
c12.dcm <- DCM(clust1, clust2, alpha = .0167)
c13.dcm <- DCM(clust1, clust3, alpha = .0167)
c23.dcm <- DCM(clust2, clust3, alpha = .0167)
c12.dcm.resid <- DCM(clust1, clust2, alpha = .0167, resid.full = TRUE)
c13.dcm.resid <- DCM(clust1, clust3, alpha = .0167, resid.full = TRUE)
c23.dcm.resid <- DCM(clust2, clust3, alpha = .0167, resid.full = TRUE)

#Sanity check: If we have random groups would we see results?
#rand.grp <- sample.int(n=3, size=ncol(vsd), replace=TRUE)
#mr1 <- ddr.vsd[,rand.grp == 1]
#mr2 <- ddr.vsd[,rand.grp == 2]
#rand.dcm <- DCM(mr1, mr2, alpha=.0167)

#######################################
######### Correlations  ############
#######################################
# Takes a DCM result as input
# Outputs a data frame 
gscors <- function(dcmres, c1 = "Correlation in C1", c2 = "Correlation in C2"){
  gs.cors <- data.frame(ClusterA_Cors = unlist(dcmres$meanCor1), 
                        ClusterB_Cors = unlist(dcmres$meanCor2))
  names(gs.cors) <- c(c1, c2)
  rownames(gs.cors) <- paste("Gene Set", 1:length(dcmres$meanCor1))
  return(gs.cors)
}
#dcm.dfs <- list(c12.dcm, c13.dcm, c23.dcm, c12.dcm.resid, c13.dcm.resid, c23.dcm.resid)
#lapply(dcm.dfs, gscors)

#######################################
######### Gene Lists  ############
#######################################
gsnames <- function(dcmres){
  dc.gene.names <- list()
  for(i in 1:length(dcmres$DC_sets)){
    dc.gene.names <- list(dc.gene.names, rownames(ddr.vsd)[dcmres$DC_sets[[1]]])
  }
  return(dc.gene.names)
}
c12.g1.genes <- rownames(ddr.vsd)[c12.dcm$DC_sets[[1]]]
c12.g2.genes <- rownames(ddr.vsd)[c12.dcm$DC_sets[[2]]]
c13.g1.genes <- rownames(ddr.vsd)[c13.dcm$DC_sets[[1]]]
c13.g2.genes <- rownames(ddr.vsd)[c13.dcm$DC_sets[[2]]]
c23.g1.genes <- rownames(ddr.vsd)[c23.dcm$DC_sets[[1]]]

write.table(paste0(c12.g1.genes, collapse=", "), file="./Results/dcm_clusters12_set1_genes.txt", 
            row.names=FALSE, col.names=FALSE, eol= "", quote = FALSE)
write.table(paste0(c12.g2.genes, collapse=", "), file="./Results/dcm_clusters12_set2_genes.txt", 
            row.names=FALSE, col.names=FALSE, eol= "", quote = FALSE)
write.table(paste0(c13.g1.genes, collapse=", "), file="./Results/dcm_clusters13_set1_genes.txt", 
            row.names=FALSE, col.names=FALSE, eol= "", quote = FALSE)
write.table(paste0(c13.g2.genes, collapse=", "), file="./Results/dcm_clusters13_set2_genes.txt", 
            row.names=FALSE, col.names=FALSE, eol= "", quote = FALSE)
write.table(paste0(c23.g1.genes, collapse=", "), file="./Results/dcm_clusters23_set1_genes.txt", 
            row.names=FALSE, col.names=FALSE, eol= "", quote = FALSE)
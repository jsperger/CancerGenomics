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
c12.dcm <- DCM(clust1, clust2, alpha = .0167)
c13.dcm <- DCM(clust1, clust3, alpha = .0167)
c23.dcm <- DCM(clust2, clust3, alpha = .0167)
c21.dcm <- DCM(clust2, clust1, alpha = .0167)
c31.dcm <- DCM(clust3, clust1, alpha = .0167)
c32.dcm <- DCM(clust3, clust2, alpha = .0167)
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
# Group 1 High Group 2 Low
c12.g1.genes <- rownames(ddr.vsd)[c12.dcm$DC_sets[[1]]]
c12.g2.genes <- rownames(ddr.vsd)[c12.dcm$DC_sets[[2]]]
# Group 2 High Group 1 Low
c21.g1.genes <- rownames(ddr.vsd)[c21.dcm$DC_sets[[1]]]
c21.g2.genes <- rownames(ddr.vsd)[c21.dcm$DC_sets[[2]]]
# Group 1 High Group 3 Low
c13.g1.genes <- rownames(ddr.vsd)[c13.dcm$DC_sets[[1]]]
c13.g2.genes <- rownames(ddr.vsd)[c13.dcm$DC_sets[[2]]]
# Group 3 High Group 1 Low
c31.g1.genes <- rownames(ddr.vsd)[c31.dcm$DC_sets[[1]]]
c31.g2.genes <- rownames(ddr.vsd)[c31.dcm$DC_sets[[2]]]
c31.g3.genes <- rownames(ddr.vsd)[c31.dcm$DC_sets[[3]]]
# Group 2 High Group 3 Low
c23.g1.genes <- rownames(ddr.vsd)[c23.dcm$DC_sets[[1]]]
# Group 3 High Group 2 Low
c32.g1.genes <- rownames(ddr.vsd)[c32.dcm$DC_sets[[1]]]
c32.g2.genes <- rownames(ddr.vsd)[c32.dcm$DC_sets[[1]]]
c32.g3.genes <- rownames(ddr.vsd)[c32.dcm$DC_sets[[1]]]
dcm.gene.list <- list(c12.g1.genes, c12.g2.genes, c21.g1.genes,
                        c21.g2.genes, c13.g1.genes, c13.g2.genes,
                        c31.g1.genes, c31.g2.genes, c31.g3.genes,
                        c23.g1.genes, c32.g1.genes, c32.g2.genes, c32.g3.genes)
gene.list.names <- c("c12.g1.genes", "c12.g2.genes", "c21.g1.genes",
                "c21.g2.genes", "c13.g1.genes", "c13.g2.genes",
                "c31.g1.genes", "c31.g2.genes", "c31.g3.genes",
                "c23.g1.genes", "c32.g1.genes", "c32.g2.genes", "c32.g3.genes")
write.genes <- function(list.of.lists, fnames){
  for(i in 1:length(fnames)){
    temp.name <- paste0("./Results/DCM/", fnames[i], ".txt") 
    write.table(paste0(list.of.lists[i], collapse=", "), file=temp.name, 
                row.names=FALSE, col.names=FALSE, eol= "", quote = FALSE)
  }
}
temp.name <- paste0("./Results/DCM/", gene.list.names[1], ".txt") 

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
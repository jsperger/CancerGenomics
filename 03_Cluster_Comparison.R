#######################################
### Nature Replicate vs DDR Clusters ##
#######################################

ddr_3_labels <- read.csv('CancerGenomics/Results/cluster_results.csv')

subtype2_labels <- read.csv('CancerGenomics/Results/subtype2_cluster_results.csv')
subtype3_labels <- read.csv('CancerGenomics/Results/subtype3_cluster_results.csv')
subtype4_labels <- read.csv('CancerGenomics/Results/subtype_cluster_results.csv')

table(subtype3_labels$NMFC3, ddr_3_labels$NMFC3)
chisq.test(table(subtype3_labels$NMFC3, ddr_3_labels$NMFC3))

table(subtype4_labels$NMFC3, ddr_3_labels$NMFC3)
chisq.test(table(subtype4_labels$NMFC3, ddr_3_labels$NMFC3))

mosaicplot(table(subtype4_labels$NMFC3, ddr_3_labels$NMFC3), shade=T, main='Comparing TCGA Subtypes and our DDR Subtypes')

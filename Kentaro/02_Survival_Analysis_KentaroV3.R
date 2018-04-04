#######################################
####### Survival Analysis ##
#######################################
# TODO: Clean data
# TODO: Figure out what time to event data to use. Both the patient and 
#       the follow up clinical info have slightly different days to death data
#       Also I think there's some progression free time data hiding somewhere
# TODO: Determine clinical variables for inclusion in the Cox model
# TODO: Get code ready

#######################################
######### Libraries ############
#######################################
library("TCGAbiolinks")
library(survival)
library("glmpath")
library("survminer")

#######################################
######### Data Preparation ######
#######################################
#Load the clinical data
clinical.query <- GDCquery(project = "TCGA-OV", 
                           data.category = "Clinical")

# TODO: Clean, figure out what's useful
clinical <- GDCprepare_clinic(clinical.query, clinical.info = "patient")
clinical.drug <- GDCprepare_clinic(clinical.query, clinical.info = "drug")
clinical.fup <- GDCprepare_clinic(clinical.query, clinical.info = "follow_up")
clinical.admin <- GDCprepare_clinic(clinical.query, clinical.info = "admin")

#Might want to only choose popular drugs
sort(table(clinical.drug$drug_name))
#INFORMATION FROMhttps://pubchem.ncbi.nlm.nih.gov/compound/paclitaxel#section=Top
# Some of the popular drugs (only 50+ uses)
# Carboplatin IS platinum based
# Taxol IS NOT
# paclitaxel IS NOT
# cisplatin  IS platinum based
# Doxil IS NOT
# Topotecan IS NOT
# Taxotere IS NOT
# Gemcitabine IS NOT
# Gemzar IS NOT

plat <- c("Carboplatin","Cisplatin")

#However, most times chemotherapy drugs are given in combos so thats the next thing to do, sort by platinum vs not use
#And Carboplatin vs Cisplatin

#Plat vs non_plat
patient_names <- names(table(clinical.drug$bcr_patient_barcode))

plat_use <- rep(NA,length(patient_names))
plat_type <- rep(NA,length(patient_names))


for (i in c(1: length(patient_names)))
{
  data <-  clinical.drug[clinical.drug$bcr_patient_barcode == patient_names[i],]
  if (any(plat %in% data$drug_name) == "TRUE"){
    plat_use[i] <- 1
  }
  else{
    plat_use[i] <- 0
    plat_type[i] <- "Neither"
  }
  if (plat[1] %in% data$drug_name & plat[2] %in% data$drug_name){
      plat_type[i] <- "Both"
  }
  else if (plat[1] %in% data$drug_name &  !(plat[2] %in% data$drug_name)){
      plat_type[i] <- plat[1]
    }
  else if (!(plat[1] %in% data$drug_name) &  plat[2] %in% data$drug_name){
      plat_type[i] <- plat[2]
    }
  else{
      plat_type[i] <- "Neither"
    }
}

temp1 <- data.frame("bcr_patient_barcode" = patient_names,plat_use,plat_type)
clinical <- merge(clinical, temp1, by="bcr_patient_barcode")

plat_use2 <- plat_use
plat_type2 <- plat_type 

plat_use <- c()
plat_type <- c()

for (i in 1:length(patient_names)){
  data <-  clinical.drug[clinical.drug$bcr_patient_barcode == patient_names[i],]
  reps <- dim(data)[1]
  if (plat_use2[i] == 1){
    plat_use <- c(plat_use, rep(1,reps))
  }
  else{
    plat_use <- c(plat_use, rep(0,reps))
  }
  
  if (plat_type2[i] == "Both"){
    plat_type <- c(plat_type, rep("Both", reps))
  }
  else if (plat_type2[i] == plat[1]){
    plat_type <- c(plat_type, rep(plat[1], reps))
  }
  else if (plat_type2[i] == plat[2]){
    plat_type <- c(plat_type, rep(plat[2], reps))
  }
  else{
    plat_type <- c(plat_type, rep("Neither", reps))
  }
}

clinical.drug <- clinical.drug
temp <- data.frame(plat_use, plat_type)
clinical.drug <- cbind(clinical.drug, temp)



#Roughly, 0.9396226 of the patients were treated with platinums. Might be enought to do a comparision study?
#Rougly 0.06 were treated with only cisplatin 
#So it seems that we have decent numbers for comparsion between platinum based and wihotut platinum studiess(al least a perliminary)



#######################################
######### Kaplan Meier ######
#######################################

#Important variables
str(clinical)
str(clinical.drug)
str(clinical.fup)



surv.vars <- c("bcr_patient_barcode", "days_to_death", "days_to_last_followup",
               "age_at_initial_pathologic_diagnosis", "race_list", "person_neoplasm_cancer_status",
               "ethnicity","neoplasm_histologic_grade", "radiation_therapy", "primary_therapy_outcome_success", "ethnicity", 
               "plat_use", "plat_type")

surv.drug.vars <- c("bcr_drug_barcode", "therapy_types", "drug_name", "regimen_indication", "plat_use", "plat_type")


surv.data <- clinical[, surv.vars]
surv.drug.data <- clinical.drug[, surv.drug.vars]
# Define the time to event. If we have the days to death, use that. 
# Otherwise use their last follow up
surv.data$EventTime <- ifelse(!is.na(surv.data$days_to_death), 
                              surv.data$days_to_death, 
                              surv.data$days_to_last_followup)
# Create an indicator variable for whether the event is death or censoring
surv.data$censored <- ifelse(!is.na(surv.data$days_to_death), 
                              FALSE, 
                              TRUE)


#Merge datasetes
cluster <- read.csv("cluster_results.csv")
cluster$NMFC3 <- as.factor(cluster$NMFC3)
surv.data <- merge(surv.data, cluster, by="bcr_patient_barcode")



# KM for drug use
km.surv.fit <- survfit(Surv(surv.data$EventTime, surv.data$censored) ~ surv.data$plat_type , conf.type = "plain")
# KM Plot for drug use
plot(km.surv.fit, main=expression(paste("Kaplan-Meier-estimate ", hat(S)[g](t), " for Different drug use")),
     xlab="t", ylab="Survival", lwd=2, col=1:4)
legend(x="topright", col=1:4, lwd=2, legend=c("Both", "Carboplatin Only", "Cisplatin Only", "Neither"))


#KM for Our clustering
km.surv.fit <- survfit(Surv(surv.data$EventTime, surv.data$censored) ~ surv.data$NMFC3 , conf.type = "plain")
plot(km.surv.fit, main=expression(paste("Kaplan-Meier-estimate ", hat(S)[g](t), " for Clusters ")),
     xlab="t", ylab="Survival", lwd=2, col=1:3)
legend(x="topright", col=1:3, lwd=2, legend=c("Cluster 1", "Cluster 2", "Cluster 3 "))
#Nearly the same survival time

#KM For Original Clustering

#######################################
######### Cox Model ######
#######################################
#Cox Models fails to find it significent.
cox.cluster <- coxph(data=  surv.data, Surv(EventTime, censored) ~ NMFC3)
summary(cox.cluster)
#Test of the proportinality assumption using Schoenfeld Residuals
test.ph <- cox.zph(cox.cluster)
test.ph
#Visual plot of the residuals over time seems like there is a significent violation of proprotional hazards
ggcoxzph(test.ph)

#Influential observations using dfbetas
ggcoxdiagnostics(cox.cluster, type = "dfbeta", var = c(NMFC3),
                 linear.predictions = FALSE, ggtheme = theme_bw())

#None of these subtypes seem paticularly more susceible to play use
cox1 <- coxph(data=  surv.data, Surv(EventTime, censored) ~ NMFC3 + plat_use + NMFC3 * plat_use)
summary(cox1)

#
cox1 <- coxph(data=  surv.data, Surv(EventTime, censored) ~ NMFC3 + plat_type )
summary(cox1)

#Are there bweteen drug differences?
cox1 <- coxph(data=  surv.data, Surv(EventTime, censored) ~ age_at_initial_pathologic_diagnosis +NMFC3 + plat_type * NMFC3)
summary(cox1)
#Type 2 seems more deadly, but if treated with CIsplatin and carboplatin it can last

#Does what drugs were used help?
cox1 <- coxph(data=  surv.data, Surv(EventTime, censored) ~ NMFC3 + plat_type * NMFC3)
summary(cox1)
#Cisplatin seems to have a very good, large effect on survival if you have subtype 2

#Is there an age thing going on?
cox2 <- coxph(data=  surv.data, Surv(EventTime, censored) ~ age_at_initial_pathologic_diagnosis +  NMFC3 + plat_type * NMFC3)
summary(cox2)
#Now subtype two seems to respond well to carboplatin and cisplatin 

#Full Model
cox2 <- coxph(data=  surv.data, Surv(EventTime, censored) ~ age_at_initial_pathologic_diagnosis + NMFC3 + ethnicity+ plat_type+ person_neoplasm_cancer_status  )
summary(cox2)
test.ph <- cox.zph(cox2)
test.ph
ggcoxzph(test.ph)
#
cox.cluster <- coxph(data=  surv.data, Surv(EventTime, censored) ~ age_at_initial_pathologic_diagnosis + NMFC3 + ethnicity + plat_type + radiation_therapy + person_neoplasm_cancer_status )
summary(cox.cluster)


# Cox model based on age. Just for example
cox.age <- coxph(data = surv.data, Surv(EventTime, censored) ~ age_at_initial_pathologic_diagnosis + race_list + neoplasm_histologic_grade + plat_type + radiation_therapy + person_neoplasm_cancer_status)
#Why am i getting nonconvergence?
summary(cox.age)
